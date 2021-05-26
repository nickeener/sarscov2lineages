#!/bin/bash

# From a reference sequence and a set of sequences belong to a specific SARS-CoV-2 lineage, identifies the variants
# that are present at high frequency in the samples and reports them as a VCF and two BED files.
# USAGE: ./lineageVariants.sh [Input Sequences] [Reference Sequence]
# Requires sarscov2phylo (https://github.com/roblanf/sarscov2phylo) environment
# Produces three files:
#	A VCF containing all variants ([file].vcf)
#	A BED containing only variants that affect amino acid sequence ([file]_aa.bed)
#	A BED containing all variants in the nucleotide form ([file]_nuc.bed)

if [ $# != 2 ]; then
    echo "usage: $0 genomes.fasta reference.fasta"
    exit 1
fi

installDir=$(dirname "${BASH_SOURCE[0]}")
input_fasta=$1
ref_fasta=$2

cores=16 # Change this value to change number of cores used for alignment step

count=$(grep -c '>' $input_fasta)
cutoff=$(echo "0.95" | bc -l) # Change this value to adjust cutoff for SNPs
cutoff_count=$(echo "$count * $cutoff" | bc -l)
cutoff_count=${cutoff_count%.*}

# Aligns sequences using sarscov2phylo alignment pipeline
$installDir/scripts/global_profile_alignment.sh -i $input_fasta -o $input_fasta.aligned -t $cores

# Converts FASTA alignment to VCF
cat $ref_fasta <(sed -re 's/^>hCoV-19\//>/' $input_fasta.aligned) \
| faToVcf stdin $input_fasta.vcf

# Removes SNPs present at less than the cutoff frequency and removes sample columns
$installDir/scripts/vcfFilter -minAc=$cutoff_count -rename $input_fasta.vcf \
| cut -f 1-8 > $input_fasta.filtered

# Makes sure there are no spaces in sample names (causes problems later if they do)
python $installDir/scripts/processFASTA.py $input_fasta

# Creates list of variants (including indels)
minimap2 --cs $ref_fasta $input_fasta.proc | paftools.js call -L 10000 - > $input_fasta.vars

# Custom script that converts SNPs from sarscov2phylo and indels from minimap2 into a single VCF and two BED files
python $installDir/scripts/lineageVariants.py -va $input_fasta.vars -vc $input_fasta.filtered -r $ref_fasta

# Remove intermediate files
rm $input_fasta.proc $input_fasta.vars $input_fasta.aligned $input_fasta.filtered $input_fasta.vcf \
    ${input_fasta}_nohead.vcf
