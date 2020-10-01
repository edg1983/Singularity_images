#!/usr/bin/env python3
'''
Script to annotate SV VCF using any tab-separated file containing var ID and a float value
Float value is annotated in INFO with the specified tag
A filter tag is added if value is < (using min) or > (using max) a specified threshold

Creator: Edoardo Giacopuzzi
18/09/2019
'''

import sys
#sys.path.append('/well/gel/HICF2/software/BRC_tools')
import os
import argparse
import gzip
import re

def parseINFO(vcfline, tag="ALLTAGS", miss="NA"):
	infos = {}
	fields = vcfline.split("\t")
	info_tags = fields[7].split(";")
	for value in info_tags:
		p = re.compile('(.+)=(.+)')
		m = p.match(value)
		if m:
			infos[m.group(1)] = m.group(2)
		else:
			infos[value] = True
	
	if tag != "ALLTAGS": 
		if tag in infos.keys():
			myinfo = infos[tag]
		else:
			myinfo = miss
	else:
		myinfo = infos
	
	return myinfo

parser = argparse.ArgumentParser(
    description='Annotate SV VCF using any tab-separated file containing var ID and a float value\n \
        If both max and min are specified they act as interval')
parser.add_argument('-v','--vcf', action='store', required=True,
                   help='VCF file with structural vars (.vcf / .vcf.gz)')
parser.add_argument('-b','--bed', action='store', required=True,
                   help='BED file for annotation')
parser.add_argument('-c','--column', action='store', required=True, type=int,
                   help='number of value column to annotate')
parser.add_argument('-i','--id', action='store', required=True, type=int,
                   help='number of var ID column')
parser.add_argument('-t','--tag', action='store', default='BEDvalue',
                   help='number of var ID column')
parser.add_argument('-l','--svlen', action='store',
                   help='maximum dimension of SV to apply the filter')
parser.add_argument('-m','--min', action='store',  type=float,
                   help='min value for filtering. If value < min variant is filtered')
parser.add_argument('-x','--max', action='store', type=float,
                   help='max value for filtering. If value > max variant is filtered')
parser.add_argument('-o','--out', action='store', required=True,
                   help='output file name')                  
args = parser.parse_args()

if args.max is None and args.min is None:
    print("At least one of --max or --min must be specified")
    exit()

filtertag = args.tag + '_filter'

bed = open(args.bed, "r")
annotation = {}

line = bed.readline()
while line:
    line = line.rstrip("\n")
    cols = line.split("\t")
    annotation[cols[args.id - 1]] = float(cols[args.column - 1])
    line = bed.readline()
bed.close()

outfile = open(args.out,"w+")

if re.search("vcf.gz$", args.vcf) is not None:
    vcf = gzip.open(args.vcf,"rt")
elif re.search("vcf$", args.vcf) is not None:
    vcf = open(args.vcf,"r")
else:
    print(".vcf or .vcf.gz extension expected for variant file")
    exit()

line = vcf.readline()
while line.startswith('##'):
    outfile.write(line)
    line = vcf.readline()

newheadlines='##INFO=<ID='+args.tag+',Number=1,Type=Float,Description="Value from '+args.bed+'">\n\
##FILTER=<ID='+filtertag+',Description="Filtered based on '+args.tag+' value">\n'
outfile.write(newheadlines)

while line.startswith('#'):
    outfile.write(line)        
    line = vcf.readline()

vcfcols = []
while line:
    filtered = False
    shortVar = True
    line = line.rstrip("\n")
    vcfcols = line.split("\t")

    svlen = int(parseINFO(line, tag="SVLEN", miss="0"))
    
    if vcfcols[2] in annotation.keys():
        value = annotation[vcfcols[2]]
        if args.svlen is not None:
            if abs(svlen) > abs(int(args.svlen)): shortVar=False 
       
        if args.max is None and args.min is not None:
            if value <= args.min and shortVar==True: filtered = True
        elif args.max is not None and args.min is None:
            if value >= args.max and shortVar==True: filtered = True
        else:
            if value < args.min or value > args.max and shortVar==True: filtered = True
    
    else:
        value = '.'

    vcfcols[7] += ';' + args.tag + '=' + str(value)
        
    if filtered == True:
        if vcfcols[6] == 'PASS' or vcfcols[6] == '.':
            vcfcols[6] = filtertag
        else:
            vcfcols[6] += ';' + filtertag
    
    outline = "\t"
    outline = outline.join(vcfcols) + "\n"
    outfile.write(outline)
    line = vcf.readline()

vcf.close()
outfile.close()
