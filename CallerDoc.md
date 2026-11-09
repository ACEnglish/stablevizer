Pipeline Description:
=====================

This is a step-by-step overview of what this pipeline does. Below are detailed descriptions of the output files.
Many details are missing and will be filled out over time and mainly in the Methods section.

### Step 1 : Find loci of any interest

We pre-select regions for analysis (at least 10bp of length variability in per-haplotype and a `scan`
interpretable plume signal or plump signal), 

### Step 2 : Build reference annotations
Intersections of the loci against SegDup, ARD, reference bed files helps identify loci which will be ignored because of
potential reference-bias

### Step 3 : Identify coverage anomalies

Donors which have unusual coverage over a locus compared to all the other loci are marked for exclusion. This results in
not every locus having 25 donors to analyze.

### Step 4 :  MeanShift cluster reads
Per-donor-haplotype, MeanShift clustering is performed to identify the germline cluster of reads (see `is_germ` below).
Note that part of this process includes the ability to 'flip' haplotypes of reads that are not `is_germ` in their Whatshap
assigned haplotype's germline distribution, but are within the other haplotype's germline distribution.

### Step 5 : Annotate Motifs
Each donor-tissue-haplotype which has at least 3 `not is_germ` reads is annotated for TR motifs using
TandemRepeatsFinder. This process starts by building a multiple-sequence alignment with abpoa of all the reads, then
building a consensus sequence using the 'is_germ' reads, then run TRF on that consensus sequence, choosing the longest
non-overlapping annotation(s) across the consensus. Now, because we have an MSA, we can map the boundaries of those
annotation(s) back to the MSA and have a way of annotating the copynumber/prurity of every read over that MSA span. This
allows us to harmonize the TR motif annotations across reads.

### Step 6 : Identify somatic outlier reads
Using the TR motif annotations, we can find reads that have ≥10bp change of the consensus sequence's motif over the
region and mark those reads as somatic.

### Step 7 : Record over-dispersion information
This signal is not fully processed (i.e. we still don't make "plump" calls), but we do record the 'flatness metrics' of
the per-donor-haplotype `is_germ` read-length distributions.

### Step 8 : Annotate instability Origin
The Origin annotator is per donor-haplotype across all samples with reference and alternate read counts identified as
those having `is_somatic` signal from the motif annotation. We also calculate a global origin from the consolidate
coverage across all donors.


Caller Results Files
====================

`meta.json`
-----------
This json file holds quick-lookup descriptions of the locus. Not all keys are guaranteed to be present depending on how
far into the processing the locus got. 

The two keys guaranteed to be present are `locus` with the `chr:start-end` key, and `ref_filters`, which is a list of
strings with possible values:
- ARD : locus overlaps an ARD region
- SegDup: locus overlaps a ucsc segdup track region
- chrX: locus is on chrX, so should be ignored until we can account for sex/ploidy

If any `ref_filters` are present, the locus is not evaluated, will have no more keys, and should be ignored.

`covanom_donors` is a list of strings with donors IDs (e.g. SMHT001) that have coverage anomalies across the cumulative
coverage from all samples. These donors are dropped from the locus and not evaluated in downstream steps.

`preclustering_read_counts` is a json seralized DataFrame of the read counts extracted by QDPI. This can be turned into
a DataFrame for easier analysis via:

```python
meta = json.load(open('meta.json'))

df = pd.DataFrame(meta['preclustering_read_counts'])
#            0    1    2  total_coverage  unphased_pct  allele_balance  highly_unphased  poor_allele_balance
#SMHT001     0  491  525            1016      0.000000        0.516732            False                False
```

These columns are:
- (index) : donor id
- 0 : number of unphased reads
- 1 : number of haplotype 1 reads
- 2 : number of haplotype 2 reads
- total_coverage : hap0 + hap1 + hap2
- unphased_pct : unphased reads / total coverage
- allele_balance : max(hap1, hap2) / (hap1 + hap2) allele balance
- highly_unphased: boolean of if unphased_pct > MAX_UNPHASED_PCT (0.50)
- poor_allele_balance : boolean of if allele_balance > MAX_ALLELE_BALANCE (0.60)

The set of donors that are highly unphased or have poor allele balance  are excluded from downstream analysis and their
ids are recorded as a list of strings in under the keys `poorphase_donors` and `poor_allelebalance_donors`

The next key is `motif_failures`, which is a list of lists containing `[donorid, protocol, haplotype]` which had
apparent length outlier signal, but failed the TRF motif annotation steps. These donor-tissue-haplotypes are presumed to
be sequencing error and/or a non-TR change causing the outlier signal. These are excluded from downstream processing and
are assumed to be non-somatic.

The next key is `unstable_protocols`, which is a count of how many times across donors a protocol was found to be
unstable. This can be turned into a pandas Series via:
```
meta = json.load(open('meta.json'))
d = pd.Series(meta['unstable_protocols'], name='count')
#3I    2
#3U    1
#3G    1
#Name: count, dtype: int64
```

The next key is `origin` which is a dictionary with keys of "protocols" and "annotation". The "protocols" holds a list
of protocols that were found to be the origin of the instability. The "annotation" value is a string. This is a 
short-cut to the `global_origin.tsv` which is a rough attempt to determine if the locus has a shared origin across all 
donors. The global origin is detailed below.

A final key of `error` may be present in some loci which failed processing. Some will have a string value of:
- "Broken reads.fasta" : The locus' reads are messed up, usually because there are full deletions over the entire locus. These loci can be ignored
- "No Plume Reads" : The locus has no outlier signal, so only the over-dispersion (plump) signal described below is
  usable
- "No Global Origin" : The experimental global origin annotation failed.


`reads.tsv`
===========
This holds all the information of the reads spanning the locus. With the exception of donors that failed the
`covanom_donors` test above. Note that some of the columns will be null depending on if the donor-tissue-haplotype made
it all the way through the pipeline

- donor : donor id
- protocol : tissue protocol code
- hap : haplotype\*. For 
- read : QDPI read name that ties to the `reads.fasta` entry name
- length : length of the read over the locus
- sample : The SMaHTid nomenclature name 
- locus : chrom:start-end of the locus this read came from
- read_name : original sequencing read name as pulled from the input bam file
- is_germ : if True, this read is within the `germline_table.tsv` boundaries of the haplotype's germline range
- germ_length : The length of the germline centroid for this donor-haplotype
- delta : germ_length - length change in this read relative to germline
- region : a comma separated list of ints for the region indexes for TR motif annotations
- motif : a comma spearated list of TR motifs, one for each region
- delta_units : Motif unit change of the motif on this read relative to germline consensus
- is_somatic : This read has a ≥10bp length deviation in the TR motif `delta_units`


Note that the columns is_germ and is_somatic are made from semi-independent processes, meaning there's no guarantee that
is_germ != is_somatic. See the pipeline description for logic on that.

Note that the `hap` _may_ be 'flipped' for reads that are outliers relative to their original haplotype assignment's
germline length, but that lie within the boundaries of the other haplotype's germline length. So hap in the `read` may
say 1, but `hap` in the table may become 2.

`germ_table.tsv`
================
This holds descriptions of the donor-tissue-haplotype germline alleles' shape as determined by the MeanShift clustering.

- donor : donor id
- hap : haplotype of reads clustered together
- nclust : number of clusters built by the MeanShift clustering
- ngerm : number of clusters assigned to the germline super-cluster
- nassign : number of reads assigned to the ngerm super-clusters
- total : total number of reads analyzed on this donor-haplotype
- bandwidth : MeanShift hyperparamter value used
- pct_germ : percent of all reads assigned to a germline super-cluster
- mean : the average length of reads assigned to germline super-cluster
- mid : median length of germline super-cluster
- lower : minimum length of germline super-cluster
- upper : maximum length of germline super-cluster
- plump_signal : Boolean of if the allele holds an overdispersion signal, detailed in pipeline description
- plume_signal : Boolean of if the allele holds an outlier signal, detailed in pipeline description

`plume_table.tsv`
=================
This holds summary information of the outlier "plume" signal from length-based calling by donor-tissue-haplotypes
Every donor-tissue-haplotype listed had MeanShift clustering signal of at least 3 reads in a single tissue that were
outliers. The table hold:

- donor : donor id
- protocol : protocol code
- hap : haplotype (1 or 2)
- False : number of non-germline reads found
- True : number of is_germline reads found

`plump_table.tsv`
=================
This holds details of the over-dispersion "plump" signals processed.
The read length distributions for each donor's cumulative coverage is analyzed for potential over-dispersion signals.
The two indicators used are 'range', where the max/min of `is_germ` read lengths is checked, and `norm_entropy`.
A donor is has its `plump_signal` set to true if:

    (range >= MIN_OVERDISPERSION_SPAN (50)) & (norm_entropy >= MIN_OVERDISPERSION_ENTROPY (0.65))

Which selects broad read-length ranges that seem reasonably flat. The full set of columns in this table are:

- donor : donor id
- hap : haplotype 
- range : max - min is_germ read length
- excess_kurtosis : measure from scipy.stats.kurtosis. Closer to 0 means more 'normal'.
- skewness : measure from stats.kurtosis.skew. More positive/negative means more skew
- p10_90_over_range : measure of p10_90_over_range over distribution of read lengths. Higher is more flat
- iqr_over_range : measure over distribution of read lengths. Higher is more flat
- norm_entropy : normalized entropy score. Higher is more flat
- plump_signal: Boolean for if there's plump signal in the distribution

`motif_annotations.tsv`
=======================
Every donor-tissue-haplotype that has an outlier signal is fed into a motif annotation procedure. See pipeline for
details. This procedure annotates the TR sequences on every read. The file `motif_annotations.tsv` holds that sequence
annotation information. Because one TR locus may have multiple adjacent tandem repeats, we can't append this information
to the `reads.tsv` because there are sometimes >1 TR annotation per-read. The columns of this file are:

- (index) : an index column that's just `range(len(rows))`.
- locus : A key of the `donor_tissue_haplotype` this read is from
- region : A number identifying the TR region of this donor_tissue_haplotype. Note two different donor-tissue-haplotypes
  may share the same region id, but aren't guaranteed to be the same TR.
- cons_start : The MSA 'germline' consensus sequence start coordinate of the TR region
- cons_end : The MSA 'germline' consensus sequence end coordinate of the TR region
- motif : TR motif identified in this region
- motif_rc : reverse compliment of the TR
- period : the length of the motif
- read : the QDPI extracted read name
- tr_len_bp : the length of the read sequence aligned within the cons_start cons_end
- copy_number : estimated copy number of the motif over the tr_len_bp of sequence
- purity : how similar the `motif * copy_number` sequence over the region is to the read sequence
- germline_alleles_bp : An estimated length of the germline allele's TR length
- nearest_allele_bp : ignore, its always equal to germline_alleles_bp
- delta_units : increase/decrease in the copy_number of this read's TR compared to the consensus 'germline'
- is_somatic : True if `abs(tr_len_bp - germline_alleles_bp) >= 10`
- is_germ : identical to the same read's `reads.tsv` `is_germ` field

`motif_table.tsv`
=================
A summary of the motif_annotations.tsv of every donor-tissue-haplotype that have `is_somatic == True` and the read counts
of germ, soma, as well as the vaf (soma / (germ + soma)


`instability_table.tsv`
=======================
This is an experimental table that attempts to summarize the `motif_annotations.tsv` somatic read instability. This is
summarized for every locus / region in the `motif_annotations.tsv` to create the columns:

- (index) : nameless range of row numbers
- locus: donor_tissue_haplotype key to the `motif_annotations.tsv` locus above
- region : same region key to above `motif_annotations.tsv`
- cons_start : same as `motif_annotations.tsv`
- cons_end : same as `motif_annotatios.tsv`
- n_reads : number of reads on the donor_tissue_haplotype used
- n_somatic : total number of `is_somatic` reads
- vaf : n_somatic / n_reads
- vaf_lb : wilson lower-bound estimate of the variant allele fraction
- frac_expansion : ToDo
- median_delta_units : median difference in TR copy_number delta units
- range_units : max minus minimum copy number units
- iqr_units : IQR of copy number units
- noise_sigma_units : ToDo
- dev_z : ToDo
- spread_z : ToDo
- burden : ToDo
- dispersion : ToDo
- instability : ToDo
- expected_shift_uni : ToDots
- protocol : ToDo
- haplotype : ToDo

`origin_table.tsv`
==================
The set of reads in the `reads.tsv` are grouped back into their original sample an counted by `is_somatic` and not
`is_somatic`. These counts represent the alternate read count and reference read count per-donor-tissue-haplotype. These
counts are then fed into the origin annotator detailed in pipeline documentation. The output of the origin annotator is:

- (index) : unnamed index that can be ignored
- idx : the "donor.haplotype" annotated
- col : the protocol_core, protocol, tissue, dedup_tissue, layer, that was determined to be the origin
- annotation : the Core/Tissue/Tissue\* Layer annotation and if it was Specific or Outlier. May be Unclassified. Also
  might be Blood
- level: in the annotation cascade that tripped to produce this annotation
- n_col : number of 'columns' in the ref/alt read counts analyzed. This is equal to the number of protocols for e.g.
  Tissue annotations or maybe 4 for the Layer level.
- n_covered : how many of the columns have (ref + alt) ≥ 10
- n_present : how many of the columns have alt ≥ 1(?)
- ONT_ref : Number of ONT reference supporting reads
- ONT_alt : Number of ONT alternate (somatic) supporting reads
- HiFi_ref : Number of HiFi reference supporting reads
- HiFi_alt : Number of HiFi alternate (somatic) supporting reads
- batch_p_value : A calculation from the batch effect test. If null, there wasn't two platforms in the origin to do the tests. if less than ~0.05, there seems to be a coverage by platform imbalance. Note that the "Outlier" annotations' batch effect test is highly suspect and needs to be improved
- vaf : variant allele fraction of the instability
- null_vaf : For "Outlier" tests, this is the null vaf modeled used to test for enrichment/depletion
- p_value : p value from the origin test
- direction : + or - for direction of the enrichment for "Outlier" tests
- reject : true if p-value < 0.05
- has_blood : True if there is alt-reads in blood (3A) protocols
- has_fibro : True if there's alt-reads in fibroblast (3AC) protocols
- n_platform : number of sequence platforms with alt-reads
- adj_pval : BH correction p-value (might not be done right because it's performed per-locus instead of across all
  origin annotations
- adj_reject : True if adj_pval < 0.05
- present_tissue_pct : percent of tissues (protocols) that have alternate read support

`global_origin.tsv`
===================
Experimental attempt to use the origin annotation infrastructure to classify the locus as a whole instead of per
donor-tissue-haplotype. This works by summing the reference/alternate reads across all samples by their
protocols-platform. This puts, for example, all Blood 3A ONT reads together to make a cumulative reference/alternate
read count. Note that when consolidating tissues, their core differences are ignored and all become "XX". So global
origin can never have a Core Specific origin.

This may be useful for quickly identifying loci that have all the alternate read support in a single tissue and
perhaps-but-just-maybe useful for finding more blatant batch effects. However, when consolidating all the read support,
most (many) of the loci are of "Outlier" origin. When populating the `meta['origin']`, only the highest VAF element is
chosen. Also, for "Outlier" results, only `direction = '+'` enrichments are considered for the `meta['origin']`

