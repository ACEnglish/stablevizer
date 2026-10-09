# stablevizer
TR instability visualization tool

## Install
```bash
python3 -m pip install .
```
Requires python ≥3.11; tested on POSIX systems; required time 2 minutes

## Usage

Example read length data is available in "example/simulated_reads.tsv"

```
stablevizer example/simulated_reads.tsv
```
Creates the output files

- `output.germline.tsv` : Per-Donor-Haplotype read length distribution and estimates of germline allele length
- `output.unstable.tsv` : Per-Donor-Tissue-Haplotype somatic unstble read length distribution
- `output.instability.png` : Instability plot of unstable tissues

Runtime appx. 1 minute for ~20k reads
