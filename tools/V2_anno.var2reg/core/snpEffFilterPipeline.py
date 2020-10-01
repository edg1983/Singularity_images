'''
@author: niko
'''
import sys
if sys.version_info[0] < 3:
    raise Exception("Requires Python 3!")

from argparse import ArgumentParser, RawDescriptionHelpFormatter, \
    ArgumentTypeError
from subprocess import *
import datetime, time
import logging
import os
import vcfpy
import re

# Necessary for including python modules from a parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.utils import *
from core.config import *

#============================================================================

pipename="SnpEff filter pipeline"


#============================================================================
usage = pipename + '''                           

  snpEff filter pipeline. 

  Copyright 2018 CCRI. All rights reserved.
  
  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0
  
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''

DEF_MAX_THREADS = 16

start_time = time.time()
success=True
    
def checksuccess(stage):
    global success
    global start_time
    elapsed_time = time.time() - start_time
    if not success:
        logging.error("Pipeline failed at stage: " + stage)
        sys.exit("Pipeline failed at stage: " +stage)
    else:
        logging.info("-----------------------------------------------------------------------------")
        logging.info("Finished stage " + stage + " in " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)) )
        logging.info("-----------------------------------------------------------------------------")
        

def runSnpEffFilter(inVcf, inFilter, outVcf, outTsv=None, filterName=None, overwrite=False ):
    config = checkTools( {"snpSift" })
    
    if outTsv is None:
        outTsv = outVcf+".tsv"
    if filterName is None:
        filterName = os.path.splitext(os.path.basename(inFilter))[0]

    start_time = time.time()
    success=True
    if overwrite or not files_exist([outVcf]):
        # ###############################################################################
        #             Filter VCF (e.g., remove non-canonical chroms, non pass variants, etc.
        # ###############################################################################
        tmp = outVcf+".tmp.vcf"
        cmd = [getTool("snpSift"), "filter"]
        cmd += ["-f", inVcf ]
        cmd += ["-e", inFilter]
        cmd += ["-i", filterName]
        cmd += ["-p"]
        success = success and pipelineStep([inVcf, inFilter], tmp, cmd, shell=True, stdout=tmp)
        if outVcf.endswith(".gz"):
            bgzip(tmp, outVcf, index=True, override=True, delinFile=True)
        else:
            os.rename(tmp, outVcf)
        checksuccess("Filtering")

        # ###############################################################################
        #             Extract fields
        # ###############################################################################
        # "ANN[*].RANK", "ANN[*].CDNA_POS", "ANN[*].CDNA_LEN", "ANN[*].CDS_POS", "ANN[*].CDS_LEN","ANN[*].AA_POS", "ANN[*].AA_LEN", "ANN[*].DISTANCE",
        headers=[]
        for h in ["CHROM", "POS", "ID", "REF", "ALT", "FILTER", "AF", "AC", "DP", "MQ",
                  "ANN[*].ALLELE", "ANN[*].EFFECT", "ANN[*].IMPACT", "ANN[*].GENE", "ANN[*].GENEID", 
                  "ANN[*].FEATURE", "ANN[*].FEATUREID", "ANN[*].BIOTYPE", "ANN[*].HGVS_C",
                  "ANN[*].HGVS_P", "ANN[*].ERRORS"]: # FIXME:  there is no INFO header for CADD_PHRED
            headers+=["\""+h+"\""]
            
        
        # add INFO fields from VCF
        reader = vcfpy.Reader.from_path(outVcf)
        pattern = re.compile("^[A-Za-z_][0-9A-Za-z_.]*\Z")
        for h in reader.header.info_ids():
            # check compatibility
            if pattern.match(h):
                headers+=["\""+h+"\""]
            else:
                logging.warn("Cannot extract INFO field " + h + " as not VCF compliant")
            
        cmd = [getTool("snpSift"), "extractFields"]
        cmd += ["-s", "\",\""]
        cmd += ["-e", "\".\""]
        cmd += [outVcf]
        cmd += headers
        success = success and pipelineStep(outVcf, outTsv, cmd, shell=True, stdout=outTsv)
        checksuccess("Extracting")



    else:
        logging.info("Skipped creation of " + outVcf + " as file already exists")
  
    return success     
        
        
##############################################################################################
#        Commandline
##############################################################################################
if __name__ == "__main__":   
    parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--inVcf", required=True, dest="inVcf", help="path to the input VCF file", type=existing_file, metavar="<VCF>")
    parser.add_argument("-f", "--inFilter", required=True, dest="inFilter", help="path to the filter txt file", type=existing_file, metavar="<TXT>")
    parser.add_argument("-c", "--config", type=existing_file, default=None, required=True, dest="config", help="Config file.")
    parser.add_argument("-o", "--outVcf", required=True, dest="outVcf", help="path to the output VCF file [default: %(default)s]", type=str, metavar="<VCF>")
    parser.add_argument("-o2", "--outTsv", required=False, dest="outTsv", help="path to the output TSV file [default: <vcf>.tsv]", type=str, metavar="<TSV>")
    parser.add_argument("-fn", "--filterName", required=False, dest="filterName", help="Filter name [default: filename of filter file]", type=str, metavar="<STR>")
    
    parser.add_argument("--overwrite", required=False, action='store_true', dest="overwrite", help="If set, the output files will be overwritten if existing")
    parser.add_argument("-l", "--log", required=False, type=str, dest="logfilename", help="optional logfile")
    
    args = parser.parse_args()
     
    # load config and merge with default tool config to enable system-specific program locations
    global config
    config = loadUserConfig(args.config)
    
    outdir=os.path.abspath(os.path.join(args.outVcf, os.pardir))
    #============================================================================
    # ensure dirs
    #============================================================================
    if not outdir.endswith("/"):
        outdir += "/"
    if not os.path.exists(outdir):
            print("Creating dir " + outdir)
            os.makedirs(outdir)
    log = os.path.join(outdir, pipename.replace(" ", "-")+".log")
    success = True    
    
    # start 
    logging.basicConfig(filename=log, level=logging.DEBUG)                
    logging.info("==========================================================")
    logging.info(pipename)
    logging.info("==========================================================")
    logging.info("Effective configuration:")
    logging.info(getConfig())
    logging.info("Started script at %s.", str(datetime.date.today()))
    

    success = success and runSnpEffFilter(args.inVcf, args.inFilter, args.outVcf, args.outTsv, args.filterName, args.overwrite)
    
    if not success:
        sys.exit("Pipeline failed")
    else:
        elapsed_time = time.time() - start_time
        logging.info("-----------------------------------------------------------------------------")
        logging.info("Finished pipeline in " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)) )
        logging.info("-----------------------------------------------------------------------------")
