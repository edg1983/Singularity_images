'''
Created on June, 2019

@author: niko
'''


from argparse import ArgumentParser, RawDescriptionHelpFormatter
import csv, sys, os

# Necessary for including python modules from a parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.utils import *

usage = '''python pipeline.py                             

  HICF2 filter pipeline

  Copyright 2019 Niko Popitsch. All rights reserved.
  
  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0
  
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''

NONE_STR = "."
VALUE_SEP = ","

DEF_gerp_cutoff = 4.0
DEF_phyloP_cutoff = 1.6
DEF_phastCons_cutoff = 0.0

IMPACTS = ["UNKNOWN", "MODIFIER", "LOW", "MODERATE", "HIGH"]
CONS = ["UNKNOWN", "LOW", "HIGH"]

#
# get ordered list of unique values
#
def getUnique(values):
    seen = set()
    seen_add = seen.add
    return [x for x in values.split(VALUE_SEP) if not (x in seen or seen_add(x))]

#
# find maximum value
#
def getMax(values, possibleValuesSorted):
    maxIdx = 0
    for v in values.split(VALUE_SEP):
        idx = possibleValuesSorted.index(v)
        if idx < 0:
            raise ValueError("Unknown value provied: " + v)
        maxIdx = max(maxIdx, idx)
    return possibleValuesSorted[maxIdx]


def toFloat(string):
    if string is None:
        return None
    if "," in string:
        string=string.split(",")[0]
    if string == ".":
        return None
    return float(string)

#  @see  http://www.academia.edu/10773044/Comparison_and_integration_of_deleteriousness_prediction_methods_for_nonsynonymous_SNVs_in_whole_exome_sequencing_studies
def getMaxCons(GERP, phyloP, phastCons):
    ccount = 0
    # gerp++ scores: negative or positive scores.
    if GERP is not None:
        if GERP > DEF_gerp_cutoff:
            ccount += 1
        else:
            ccount -= 1

    # PhyloP score: negative or positive scores, range about -7.. +7.
    if phyloP is not None:
        if phyloP > DEF_phyloP_cutoff:
            ccount += 1
        else:
            ccount -= 1

    # phastCons: [0, 1]: probability that position is undernegative selection (i.e., conserved)
    if phastCons is not None:
        if phastCons > DEF_phastCons_cutoff:
            ccount += 1
        else:
            ccount -= 1
    if ccount > 0:
        return "HIGH"
    elif ccount < 0:
        return "LOW"
    else:
        return "UNKNOWN"


def getWGS500AF(TC, HM, HT):
    if TC is None or HM is None or HT is None:
        return None
    AC = 2 * float(HM.split(",")[0]) + float(HT.split(",")[0])
    AN = 2 * float(TC.split(",")[0])
    if AN == 0:
        return None
    return AC / AN

RARENESS = ["UNKNOWN", "VERY_RARE", "RARE", "COMMON", "VERY_COMMON"]
RARENESS_CUTOFF = [ 0, 0.0001, 0.001, 0.02, 1] 

def getMaxPopAF(A1000G_AF, UK10K_AF, gnomAD_AF, ExAC_AF, WGS500_AF):
    af = 0
    m = []
    if A1000G_AF is not None:
        if A1000G_AF > af:
            af = A1000G_AF
        m += ["1000G=" + str(A1000G_AF)]
    if UK10K_AF is not None:
        if UK10K_AF > af:
            af = UK10K_AF
        m += ["UK10K=" + str(UK10K_AF)]
    if gnomAD_AF is not None:
        if gnomAD_AF > af:
            af = gnomAD_AF
        m += ["gnomAD=" + str(gnomAD_AF)]
    if ExAC_AF is not None:
        if ExAC_AF > af:
            af = ExAC_AF
        m += ["ExAC=" + str(ExAC_AF)]
    if WGS500_AF is not None:
        if WGS500_AF > af:
            af = WGS500_AF
        m += ["WGS500=" + str(WGS500_AF)]
        
    ret = []
    for idx, val in enumerate(RARENESS_CUTOFF):
        if af <= val:
            ret = [RARENESS[idx], ",".join(m)]
            break
        
    return ret
    

ROIS = ["RG", "LowComplexity", "LowMappability", "TopVariableGenes", "EncodeReg", "PanelApp"]
FIELDS = ["CHROM", "POS", "REF", "ALT", "FILTER", "ID", "Genes", "Effects", "CADD_PHRED", "maxImpact", "maxPopAF", "maxPopAF2", "maxCons", "GeneIds", "ROI"]

def postfilter(infile, outfile, overwrite=False):
    if files_exist(outfile):
        if overwrite:
            removeFile(outfile)
        else:
            print("Outputfile exists!")
            return False
    with open(outfile, 'a') as out:
        # write header
        out.write("\t".join(FIELDS) + "\n")
        with open(infile, 'rt') as tsvin:
            tsvin = csv.reader(tsvin, delimiter='\t')
            headers = []
            for row in tsvin:
                if (row[0].startswith("CHROM")):
                    headers = row
                    continue
                if (row[0].startswith("#")):
                    continue
                d = {}
                for i in range(len(row)):
                    d[headers[i]] = row[i]
                o = {}
                for f in FIELDS:
                    if f == "Effects":
                        o[f] = getUnique(d.get("ANN[*].EFFECT"))
                    elif f == "Genes":
                        o[f] = getUnique(d.get("ANN[*].GENE"))
                    elif f == "GeneIds":
                        o[f] = getUnique(d.get("ANN[*].GENEID"))
                    elif f == "maxImpact":
                        o[f] = getMax(d.get("ANN[*].IMPACT"), IMPACTS)
                    elif f == "maxCons":
                        o[f] = getMaxCons(toFloat(d.get("GERP")), toFloat(d.get("phyloP")), toFloat(d.get("phastCons")))
                    elif f == "ROI":
                        rois = []
                        for r in ROIS:
                            if d.get(r) is not None and d.get(r).lower() == "true":
                                rois += [r]
                        if not rois:
                            rois=["."]
                        o[f] = ",".join(rois)
                    elif f == "maxPopAF":
                        continue
                    elif f == "maxPopAF2":
                        popAF = getMaxPopAF(toFloat(d.get("A1000G_AF")),
                                        toFloat(d.get("UK10K_AF")),
                                        toFloat(d.get("gnomAD_AF")),
                                        toFloat(d.get("ExAC_AF")),
                                        getWGS500AF(d.get("WGS500_TC"), d.get("WGS500_HM"), d.get("WGS500_HT")))
                        o["maxPopAF"] = popAF[0]
                        o["maxPopAF2"] = popAF[1]
                    else:
                        o[f] = d.get(f)
                
                
                first = True
                for f in FIELDS:
                    if not first:
                        out.write("\t")
                    first = False
                    if (type(o[f]) is list):
                        out.write(VALUE_SEP.join(o[f]))
                    else:
                        out.write(str(o[f]))
                out.write("\n")
    print("Finished.")           
    return True
                
##############################################################################################
#        Commandline
##############################################################################################
if __name__ == "__main__":   
    #============================================================================
    parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--in", type=existing_file, required=True, dest="infile", metavar="infile", help="Input TSV file")
    parser.add_argument("-o", "--out", type=str, required=True, dest="outfile", help="output TSV file")
    parser.add_argument("--overwrite", required=False, action='store_true', dest="overwrite", help="If set, the output files will be overwritten if existing")
    args = parser.parse_args() 
    #============================================================================
    postfilter(args.infile, args.outfile, args.overwrite)
