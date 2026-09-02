# pastri
Pattern Analysis for Somatic Tandem Repeat Instability

Install via
```bash
python3 -m pip install .
```

Usage See `pastri -h` for available commands

# Data Setup
Currently, pastri assumes that all metadata can be tracked through SMaHT nomenclature sample ids. Future release will
allow definitions of custom donor-tissue-platform-etc fields.

Use `qdpi` to extract read deltas, and optionally sequences, from long read alignments over tandem repeat regions.
```bash
Link to qdpi and show an example command
```
Give full paths to the sets of files

```bash
ls `pwd`/*.qdpi.bed.gz > qdpi_files.txt
ls `pwd`/*.fa > fasta_files.txt
```

Then, use `pastri tr` to pull relevant read information from the files for a particular locus
```bash
pastri tr -o output.reads.tsv chrom:start-end -q qdpi_files.txt -f fasta_files.txt
```

# Plume
When a central mass of read lengths are clustered around a germline TR length, any remaining
outlier read lengths are deemed somatic 'plumes'. This pattern can be searched for in a set of reads over a locus with:
```bash
pastri plume output.reads.tsv -o plume ...
```

# Origin
Origin annotation takes the plume `output.anno_reads.tsv` and searches for the simplest core/tissue/layer origin of the
non-germline reads. 

```bash
pastri origin -o origin.tsv plume.anno_reads.tsv
```
Documentation due

# MSA
Documentation due

# Motif
Documentation due

# KSOM

Kmer-based Self-Organizing Map.

### Usage:

1. Start by building the kmer vectors from reference sequences within a bed file
```bash
pastri ksom build-kvec -b input.bed -f reference.fa -o kvec.npz
```

2. Build a SOM from those kmer vectors with
```bash
pastri ksom build-som -o som.pkl -p kvec.npz
```

3. Map the regions to the SOM via
```bash
pastri ksom map-bed -b input.bed -s som.pkl -o out.map.bed
```
The `out.map.bed` will hold each input bed line's first 3 columns along with two additional `X` and `Y` columns.

It is also possible to directly map kvec files to the SOM via:
```bash
pastri ksom map-kvec -k kvec.npz -s som.pkl -o out.map
```
However, this doesn't track bed-like locus information, so other infrastructure tracking the kvec indices would be
needed

4. Plot the SOM
```bash
pastri ksom plot-som -b out.map.bed -s som.pkl -o general_count.png
```

5. Generate sequence stats and plot
Coloring the SOM by e.g. (G+C)% is often more informative. We can generate those stats via:
```bash
pastri ksom seqstat-bed -b input.bed -r reference.fa -o seqstat.bed
```

We can then combine those stats with the map to make it easier to push the relevant columns to plot-som
```bash
paste out.map.bed <(cut -f4- seqstate.bed) > out.seqstat.map.bed
```

And then color the SOM plot by mean (G+C)% via:
```bash
pastri ksom plot-som -b <(cut -f1-6 out.seqstat.map.bed) \
    -s som.pkl \
    -o gcpct.png \
    --metric mean \
    -T "(G+C)%", --title "GC Percent"
```

Or color the SOM by homopolymer percent via:
```bash
pastri ksom plot-som -b <(cut -f1-5,7 out.seqstat.map.bed) \
    -s som.pkl \
    -o hompct.png \
    --metric mean \
    -T "Homopolymer%" --title "TR Region Homopolymers"
```

6. Subset Enrichment Tests
Test if a subset of TR regions is significantly enriched in a subset of SOM neurons.

```bash
pastri ksom enrichment -s .som.pkl -n out.map.bed -m subset.map.bed -o subset.test.tsv
```

Then, visualize the test results against the SOM with
```bash
pastri ksom plot-som -b <(cut -f1-5,7 out.seqstat.map.bed) \
    -s som.pkl \
    -o Test.png \
    --enrichment subset.test.tsv \
    --metric mean \
    --bar-title "Homopolymer%" --title "Enrichment Test" \
```
