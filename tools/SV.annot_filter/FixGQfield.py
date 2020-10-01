#!/usr/bin/env python3
from collections import OrderedDict
import argparse
import gzip
import sys, os

def checkFile(file):
    if not os.path.isfile(file):
        sys.exit("CRITICAL! " + file + " do not exists!")
    
def parseFORMAT(format_string, sample_string):
    vcfFORMAT = format_string.split(":")
    vcfSAMPLE = sample_string.split(":")
    values = OrderedDict(zip(vcfFORMAT, vcfSAMPLE))
    return values

def updateFORMAT(sample_values, field, old_string, new_string):
    updated = 0
    if sample_values['GT'] != './.':
        if sample_values[field] == old_string:
            updated = 1
            sample_values[field] = new_string
    return updated, sample_values

def roundFORMAT(sample_values, field):
    updated = 0
    if field in sample_values.keys() and sample_values[field] != '.':
        updated = 1
        sample_values[field] = round(float(sample_values[field]))
    return updated, sample_values

def tokenize(line,sep):
    line = line.rstrip('\n')
    line = line.split(sep)
    return line

parser = argparse.ArgumentParser(
    description='Round GQ values to make them integers, so they can be processed dy bcftools')
parser.add_argument('-v','--vcf', action='store', required=True,
                   help='VCF file with structural vars (.vcf / .vcf.gz)')
parser.add_argument('-o','--out', action='store', required=True,
                   help='output file name')                  
args = parser.parse_args()

checkFile(args.vcf)
if args.vcf.endswith("vcf.gz"):
    vcf = gzip.open(args.vcf,"rt")
elif args.vcf.endswith("vcf"):
    vcf = open(args.vcf,"r")

outfile = open(args.out, "w+")
field = 'GQ'

line = vcf.readline()
while line.startswith('#'):
    outfile.write(line)
    line = vcf.readline()

nvars = 0
ngenos = 0
corrected = 0
while line:
    nvars += 1
    line = tokenize(line, "\t")
    format_string = line[8]
    for i in range(9, len(line)):
        ngenos += 1
        sample_values = parseFORMAT(format_string,line[i])
        is_updated, new_values = roundFORMAT(sample_values,field)

        #is_updated, new_values = updateFORMAT(sample_values,field,old_value,new_value)

        line[i] = ":".join([str(x) for x in new_values.values()])
        corrected += is_updated
        
    outfile.write("\t".join(line) + "\n")
    line = vcf.readline()

outfile.close()
vcf.close()   
print(nvars, " variants")
print(ngenos, " genotypes")
print(corrected, " fixed ", field, " values")
