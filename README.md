# Identifying Complete Constellation of Variants in SARS-CoV-2 Lineages

## We Need All the Variants

The massive number of SARS-CoV-2 genomes available on public databases such as GISAID make it possible for us to track the genetic lineage of the virus, leading to several major viral lineages to be identified (B.1.1.7, B.1.351, etc.). Very often, these lineages are identified based on a small set of characteristic variants that are present in nearly every sample from that lineage and are of biological importance (vaccine escape, affects transmissibility, ect.). But this is ignoring a significant chunk of the true variation in these viral lineages as each also has many fixed variants that have either unknown relevance or are known to be inconsequential. Knowing this full set of variants is vital as when new lineages crop up that have similar sets of variants with known lineages, we may only be able to tell them apart due to their differences in variants that are not part of the common characteristic set.

## Lineage Variants 

This repo contains a set of scripts for identifying all high frequency variants in a set of SARS-CoV-2 genomes and converting these to BED files that can be used to create a custom track in the SARS-CoV-2 UCSC genome browser (https://genome.ucsc.edu). It will also create a VCF in standard format describing the variation.

### Usage

Clone this repository using: ```git clone https://github.com/nickeener/sarscov2lineages.git```

To ensure you have all the required dependencies, a .yml file for creating a custom conda environment with all the correct dependencies has been provided (this requires you to already have some form of anaconda installed, see: https://www.anaconda.com/). Simply call ```conda env create -f environment.yml``` to create the environment and then ```conda activate sarscov2lineages``` to activate the environment.

Depending are your system, you can change the number of cores you wish to use for the multiple alignment step of the pipeline by editing line 12 of lineageVariants.sh to your desired number of cores (more cores means it will run faster).

Once you have the necessary dependencies and have set your optimal number of cores the pipeline can be run using the following command:

```./lineageVariants.sh [path to multi-FASTA containing assembled genomes] [path to SARS-CoV-2 reference FASTA]```

(The SARS-CoV-2 reference is also included in this repo as the file NC_045512.fa)

