#!/bin/bash

# From a reference sequence and a set of sequences belong to a specific SARS-CoV-2 lineage, identifies the variants
# that are present at high frequency in the samples and reports them as a VCF and two BED files.
# USAGE: ./lineageVariants.sh [Input Sequences] [Reference Sequence]
# Requires sarscov2phylo (https://github.com/roblanf/sarscov2phylo) environment
# Produces three files:
#	A VCF containing all variants ([file].vcf)
#	A BED containing only variants that affect amino acid sequence ([file]_aa.bed)
#	A BED containing all variants in the nucleotide form ([file]_nuc.bed)

cores=16 # Change this value to change number of cores used for alignment step

count=$(grep -c '>' $1)
cutoff=$(echo "0.95" | bc -l) # Change this value to adjust cutoff for SNPs
cutoff_count=$(echo "$count * $cutoff" | bc -l)
cutoff_count=${cutoff_count%.*}

# Aligns sequences using sarscov2phylo alignment pipeline
scripts/global_profile_alignment.sh -i $1 -o $1.aligned -t $cores

# Converts FASTA alignment to VCF
cat $2 <(sed -re 's/^>hCoV-19\//>/' $1.aligned) | scripts/faToVcf stdin $1.vcf

# Removes SNPs present at less than the cutoff frequency and removes sample columns
scripts/vcfFilter -minAc=$cutoff_count -rename $1.vcf | cut -f 1-8 > $1.filtered

# Makes sure there are no spaces in sample names (causes problems later if they do)
python scripts/processFASTA.py $1

# Create list of variants (including indels)
minimap2 --cs $2 $1.proc | paftools.js call -L 10000 - > $1.vars

# Custom script that converts SNPs from sarscov2phylo and indels from minimap2 into a single VCF and two BED files
python scripts/lineageVariants.py -va $1.vars -vc $1.filtered -r $2

# Remove intermediate files
rm $1.proc $1.vars $1.aligned $1.filtered $1.vcf $1_nohead.vcf

