# Usage

(Optional) Generate sequence information of the regions and do a greedy furthest-point downsampling. This is an attempt
to even out the sequence-context representations before building SOM in order to avoid overly dense clusters. Unknown
how good it works at evening out, but it does speed up SOM building.

```bash
python ksomtr.py bed-seqstat \
    -b models/adotto.v2.trgt.lite.bed \
    -r ~/code/references/grch38/GRCh38_1kg_mainchrs.fa \
    -o models/adotto.v2.trgt.lite_seqstat.bed

python ksomtr.py downsample -n 100000 \
    -b models/adotto.v2.trgt.lite_seqstat.bed \
    -o models/adotto.v2.trgt.lite_seqstat_down_sampled.bed 
```

Okay, maybe don't do this. I noticed that the downsampled is uniform across the SOM, but when I apply, somehow there's a
couple of neurons that have a whole lot of the points. So I need to train on the full if I want a flat map.

Vectorize reference sequences
```bash
python ksomtr.py kvec --normalize \
    -b models/adotto.v2.trgt.lite_seqstat_down_sampled.bed \
    -r ~/code/references/grch38/GRCh38_1kg_mainchrs.fa \
    -o models/downsample.normalize.kvec.npz
```

Build the SOM
```bash
python ksomtr.py som -p -i 400000 models/downsample.normalize.kvec.npz -o models/downsample.normalize.som.pkl
```

Map the regions to the SOM to get their X/Y coords
```bash
python ksomtr.py bed-map -s models/downsample.normalize.som.pkl \
    -r ~/code/references/grch38/GRCh38_1kg_mainchrs.fa \
    -b models/adotto.v2.trgt.lite.bed > adotto.full_mapping.bed
```

### TODO ###

- Allow sequences (-f) instead of (-b/-r)
  - And abstract out that argparse adding
- Do I need to do so many iter?
- Save optimization in the som output

Plot how many TR loci are in each point of the SOM
```bash
python ksomtr.py plot -s models/downsample.normalize.som.pkl \
    -b adotto.full_mapping.bed \
    -t "TR Count" \
    -o plot.png
```

Plot where the pathogenic TR occur
```bash
bedtools intersect -u -a adotto.full_mapping.bed -b pathogenic.bed > pathogenic_mapping.bed

python ksomtr.py plot -s models/downsample.normalize.som.pkl \
    -b adotto.full_mapping.bed \
    -t "Pathogenic" \
    -o patho_plot.png
```

- Do I want to improve the auto plotting? Or even do it at all?
- I still want to be able to put in markers
- I really wish the bed files 'carried' information forward - could move region_parse to pandas
- Expose every method. I mean, why not be able to calc hompolymer from command line


# Notes:

I need to do a better job of putting in sites
Like, Instead of putting in all TRs, what If I sample some kind of even subset?
I want the map to be as flat as possible and have good spread.
So how do we pick representative kmers...?

### Kmer Self-Organizing Map for TR Viz

Okay - I got these little tools
I don't need to pipeline all of this. 
I just need to clean up the plotting code so that it's easier to run more generally.
So go get familiar with that thing, I reckon.

1. Command that takes a bed file and double-strand kmer-vectorizes the reference sequences spanned
    - Tool (notebook) that will build the map from that
      - This should be inside the command. So the output is going to be .. joblib I guess of
```json
      { 
        'som' : minisom object that allows mapping of new sequences
        'gc_pct': x/y coordinates of the som to an average GC percent. Will be default coloring
        'locus_lookup': locus mapping
        'region_lookup': intersection lookup
    }
```
    Is that sufficient for doing all the plotting you want to do? 
    Like, quickly mapping plots is a main goal of this thing.

    Nice to have - in addition to gc_pct, you could do patho percent or homopolymer percent?
    Many of those should be custom lookups, I think, perhaps don't save too much

2. So, what commands do you want?

    ksomtr build bed_file reference --mini --som --parameters?

    ksomtr bed-map bed_file som_file --plot 
        Outputs the plot of just heatmap
        Outputs the annotated bed_file of the X/Y coordinates it goes to.

    ksomtr bed-map bed_file som_file --plot --mark
        Adds marks instead of / in addition to heatmap?
        I guess a default heatmap of GCpct? So record that

    ksomtr seq-map given 
        given sequences, do a mapping into the som - so we build the kmervectors
        
    And keep the plotting minimal, but expose from library easy hooks into it for custom plotting

Tool that will take another bed file (or a set of locus ids) and will heatmap the counts
