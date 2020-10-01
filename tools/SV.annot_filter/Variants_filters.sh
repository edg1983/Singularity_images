#!/bin/bash

# Filters for SVs as suggested in
# https://www.biorxiv.org/content/10.1101/508515v1
# This is set up to run using singularity images containing some supporting files in /resources

input=$1 #input file name without path
output=$2 #output file name without path
workingdir=$3 #Where input file is present. Will be bound to singularity /data folder
SVpipeline_img=$3 #Path to SVpipeline.sif
SVannot_img=$4 #Path to SV.annot_filter.sif

excluderegions="/resources/hg38_CNV_ExcludeRegions+LowComplexity.bed"
svtools="singularity exec --bind $workingdir:/data $SVpipeline_img svtools"
SVannot="singularity exec --bind $workingdir:/data $SVannot_img"

# First clean the file removing vars with incomplete stats due to high missingness
echo "1. Remove missing calls"
(zgrep "#" $workingdir/$input && zgrep -v "#" $workingdir/$input | awk '$9 != "GT" && $9 != "GT:CN" {print; }') > $workingdir/temp_nomiss.vcf

# Fix the GQ values
echo "2. Fix GQ field"
$SVannot FixGQfield.py -v /data/temp_nomiss.vcf -o /data/temp_fixGQ.vcf

echo "3. Apply filters"
# 1. For small deletions < 1000bp
# LowConfSmallDEL filter: the variant doesn't have split-read support in any sample (AS == 0 in all samples) 
$SVannot bcftools filter -m + -s LowConfSmallDEL -e 'INFO/SVTYPE == "DEL" && INFO/SVLEN >= -1000 && N_PASS(AS > 0) == 0' /data/temp_fixGQ.vcf > $workingdir/temp_1.vcf

# 2. For INV
# INV_LowMSQ: MSQ < 150
# INV_LowLumpyEvidence: neither split-read nor paired-end lumpy evidences are at least 10% of total evidence 
$SVannot bcftools filter -m + -s INV_LowMSQ -e 'INFO/SVTYPE == "INV" && INFO/MSQ < 150' /data/temp_1.vcf > $workingdir/temp_2.vcf
$SVannot bcftools filter -m + -s INV_LowLumpyEvidence -e 'INFO/SVTYPE == "INV" && ((INFO/SR / INFO/SU) < 0.1 || (INFO/PE / INFO/SU) < 0.1)' /data/temp_2.vcf > $workingdir/temp_3.vcf

# 3. For BND
# BND_LowMSQ: MSQ < 250 
$SVannot bcftools filter -m + -s BND_LowMSQ -e 'INFO/SVTYPE == "BND" && INFO/MSQ < 250' /data/temp_3.vcf > $workingdir/temp_4.vcf

# 4. Size filtering to remove DEL/DUP < 50 bp (better captured by small-variant calling)
# VerySmallSV
$SVannot bcftools filter -m + -s VerySmallSV -e '(INFO/SVTYPE == "DEL" && INFO/SVLEN > -50) || (INFO/SVTYPE == "DUP" && INFO/SVLEN < 50)' /data/temp_4.vcf > $workingdir/temp_5.vcf

# 5. Filter variants that have multiple overlap with LowComplexity and ExcludeRegions
# SMALL SV (< 1kb) with overlap >70% are filtered
# First convert SV VCF to BED (BND variants are discarded)
$svtools vcftobedpe -i /data/temp_5.vcf -o /data/temp_5.bedpe
$svtools bedpetobed12 -i /data/temp_5.bedpe -o /data/temp_5.bed
grep -v "BND" $workingdir/temp_5.bed | tail -n+2 | cut -f1-4 | tr ";" "\t" | cut -f1-5 | sed 's/ID=//g' > $workingdir/temp_5.noBND.forOverlap.bed
#Compute the overlap with RLCRs regions
$SVannot bedtools coverage -wa -a /data/temp_5.noBND.forOverlap.bed -b $excluderegions > $workingdir/temp_5.noBND.overlapExcludeRegions.txt
#Apply filter based on 70% threshold (using variant ID as reference)
$SVannot bgzip /data/temp_5.vcf
$SVannot FilterByValue.py -v /data/temp_5.vcf.gz -b /data/temp_5.noBND.overlapExcludeRegions.txt -c 9 -i 5 -t LowComplexity -x 0.7 -l 1000 -o /data/$output

#REMOVE TEMP FILES
echo "4. Remove temp files"
rm $workingdir/temp_*

echo "### ALL DONE! ###"
echp "Output file is: $workingdir/$output"

# 6. Additional cohort analysis
# Remove outlier samples
# samples where variant counts (for any SV type) were >6 median absolute deviations from the median count for that type
