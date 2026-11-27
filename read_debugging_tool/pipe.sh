
# Parameters
# Set the refeference
ref=/Users/English/code/references/grch38/GRCh38_1kg_mainchrs.fa

# Set the region of interest
chrom=chr1
start_pos=99682978
end_pos=99683197
# How much buffer upstream/downstream to build the flanks
BUFF=200

# prefix to output files
name=$1
# Aligned long-read BAM
bam=$2
SRCDIR=$( dirname "${BASH_SOURCE[0]}" )

up_pos=$(expr $start_pos - $BUFF)
dn_pos=$(expr $end_pos + $BUFF)

echo -e "${chrom}\t${start_pos}\t${end_pos}" > region.bed
samtools view --region-file region.bed -O BAM ${bam} | samtools fasta > reads.fa
samtools view --region-file region.bed -O BAM ${bam} | samtools fastq > reads.fq

# Extract Flanks
samtools faidx $ref ${chrom}:${up_pos}-${start_pos} | sed 's/>.*/>region1/' > up_flank.fa
samtools faidx $ref ${chrom}:${end_pos}-${dn_pos} | sed 's/>.*/>region1/' > dn_flank.fa

# Manually rename the flank fasta regions to the same name
minimap2 -x map-ont -c -k11 -w3 -n2 -m25 -s25 -p 0.5 -N 50 up_flank.fa reads.fa > up_flank.paf
minimap2 -x map-ont -c -k11 -w3 -n2 -m25 -s25 -p 0.5 -N 50 dn_flank.fa reads.fa > dn_flank.paf

python ${SRCDIR}/extract_between.py reads.fa up_flank.paf dn_flank.paf --inserts ${name}.inserts.fa --clipped ${name}.clipped.fa
python ${SRCDIR}/plot_read_quality.py -i ${name}.inserts.fa -c ${name}.clipped.fa -p ${name} reads.fq
