#!/bin/bash
set -euxo pipefail

## needs GATK to be installed in your path

# intersect high confidence regions with target selection bed file of the kit
bedtools intersect -a no_chr_highconf.bed -b Agilent_v7.bed > better_overlap.bed

# indexing
gatk IndexFeatureFile -I no_chr_ref.vcf
gatk IndexFeatureFile -I no_chr_query.vcf

# some filtering
gatk SelectVariants -R ../genome.fasta --exclude-non-variants --remove-unused-alternates --intervals better_overlap.bed -V no_chr_ref.vcf -O clean_no_chr_ref.vcf
gatk SelectVariants -R ../genome.fasta --exclude-filtered --exclude-non-variants --remove-unused-alternates --intervals better_overlap.bed -V no_chr_query.vcf -O clean_no_chr_query.vcf

# left-alignment & trimming
gatk LeftAlignAndTrimVariants -R ../genome.fasta -V clean_no_chr_query.vcf -O final_clean_no_chr_query.vcf

# calculate concordance
gatk Concordance -R ../genome.fasta --eval final_clean_no_chr_query.vcf --truth clean_no_chr_ref.vcf -S summary.tsv --intervals better_overlap.bed
