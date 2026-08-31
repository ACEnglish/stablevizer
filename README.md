# pastri
TR instability visualization tool

Documentation pending.

Install via
```bash
python3 -m pip install .
```

Usage See `stablevizer -h`


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

3. Map the kvec to the SOM via
```bash
pastri ksom map-kvec -k kvec.npz -s som.pkl -o out.map
```

It is also possible to directly map bed files to the SOM via:
```bash
pastri ksom map-bed -b input.bed -s som.pkl -o out.map.bed
```
The `out.map.bed` will hold each input bed line's first 3 columns along with two additional `X` and `Y` columns.

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
