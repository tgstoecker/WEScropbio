#!/bin/bash
set -euxo pipefail

## tis script adds a GT field for a varlociraptor output

# left-align
bcftools norm -m-any targets_filt_norm_final.vcf > la_norm_final.vcf

# vembrane table GT field
vembrane table --header 'CHROM,POS,GT' 'CHROM, POS, min( [ (INFO["PROB_HET"], "0/1"), (INFO["PROB_HOM"],"1/1"), (INFO["PROB_ABSENT"], "0/0") ], key = lambda x: x[0])[1] ' la_norm_final.vcf > la_norm_final_GT.table

# reduce to just the GT column without header
tail -n +2 la_norm_final_GT.table | awk -F "\t" '{print $3}' > la_norm_final_just_GT

# GT header line - add manually
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">

#adding GT to FORMAT field
#grep -v '^#' la | awk -F "\t" '{print "GT:"$9}'
grep -v '^#' la_norm_final.vcf | awk -F "\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,"GT:"$9,$10}' OFS="\t" > la_norm_final_GT_field_no_header

# get just the header
grep '#' la_norm_final.vcf > la_norm_final_just_header

# merge GT column into table
paste -d "\t" la_norm_final_GT_field_no_header la_norm_final_just_GT | sed 's/\r$//' | awk -F "\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11":"$10}' OFS="\t" > la_norm_final_GT_field_values_no_header

# concatenate header and merged table
cat la_norm_final_just_header la_norm_final_GT_field_values_no_header > la_norm_final_GT_field_values_with_header


# add GT FORMAT header line:
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
awk '!found && /FORMAT/ { print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"; found=1 } 1' la_norm_final_GT_field_values_with_header > la_norm_final_GT_field_values_with_header.vcf
