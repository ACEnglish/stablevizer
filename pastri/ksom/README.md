# Usage

(Optional) Generate sequence information of the regions and do a greedy furthest-point downsampling. This is an attempt
to even out the sequence-context representations before building SOM in order to avoid overly dense clusters. Unknown
how good it works at evening out, but it does speed up SOM building.

Okay, maybe don't do this. I noticed that the downsampled is uniform across the SOM, but when I map the full catalog, a
couple neurons have way too many points and the TR density isn't very flat. Maybe something interesting in that
happening w.r.t. certain sequence contexts being over represented in the TR space? IDK, but I moved to full catalog SOM


```bash
python ksomtr.py bed-seqstat \
    -b models/adotto.v2.trgt.lite.bed \
    -r ~/code/references/grch38/GRCh38_1kg_mainchrs.fa \
    -o models/adotto.v2.trgt.lite_seqstat.bed

python ksomtr.py downsample -n 100000 \
    -b models/adotto.v2.trgt.lite_seqstat.bed \
    -o models/adotto.v2.trgt.lite_seqstat_down_sampled.bed 
```


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
Slightly faster option where you don't regenerate the kvec
```bash
python ksomtr.py kvec-map -s models/downsample.normalize.som.pkl \
    -k models/downsample.normalize.kvec.npz \
    -o models/downsample.normalize.map.txt

paste models/adotto.v2.trgt.lite_seqstat_down_sampled.bed models/downsample.normalize.map.txt > adotto.full_mapping.bed
```
Then you can `sort | uniq` the map to get an idea of the spread of loci. Or go to the plotting Notebook


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
### TODO ###

- Allow sequences (-f) instead of (-b/-r)
  - And abstract out that argparse adding
- Do I need to do so many iter?
- I really wish the bed files 'carried' information forward - could move region_parse to pandas
- I think I can merge seq_to_kvec and seq_to_kmer now that it is out of kanpig
- Do I want to improve the auto plotting? Or even do it at all?
- I still want to be able to put in markers
- Make a `build` command that does all the above steps automatically 
- Expose every method. I mean, why not be able to calc hompolymer from command line
- Revisit entropy - maybe pull from truvari.
- Better structure? Think plink, the files sit in a directory so they can talk to each other auto-magically
- Packing - could go beside e.g. `stablevizer inst` as `stablevizer som` But then you got commands like `stablevizer som
  som`, which is gross
