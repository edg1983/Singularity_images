'''
@author: niko
'''

from argparse import ArgumentParser, RawDescriptionHelpFormatter, \
    ArgumentTypeError
import datetime, time
import logging
from subprocess import *
import sys, os, json

# Necessary for including python modules from a parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.utils import *

usage = '''python pipeline.py                             

  vcfanno pipeline

  Copyright 2019 Niko Popitsch. All rights reserved.
  
  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0
  
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''
DEF_MAX_THREADS = 4


def runVcfanno(inVCF, outVCF, configF, luaF=None, threads=DEF_MAX_THREADS, overwrite=True, additionalParameters=[], settings=[]):
    
    checkTools({"vcfanno" })
    
    success = True
    
    if files_exist(outVCF) and not overwrite:
        print("VCF file " + inVCF + " already exists! SKIPPING re-creation...")
        logging.warn("VCF file " + inVCF + " already exists! SKIPPING re-creation...")
        return success
    
    iscompressed = outVCF.endswith('.gz')
    vcfraw = outVCF + ".raw.vcf"

    cmd = settings + [ getTool("vcfanno") ]
    if luaF is not None:
        cmd += ["-lua", luaF ]
    cmd += [ "-p", str(threads) ]
    cmd += [ configF ]
    cmd += [ inVCF ]
    success = success and pipelineStep(inVCF, vcfraw, cmd, shell=True, stdout=vcfraw)    
    
    if iscompressed:
        bgzip(vcfraw, outFile=outVCF, index=True, delinFile=True)
    else:
        os.rename(vcfraw, outVCF)
    
    return success


if __name__ == "__main__":
    # ============================================================================
    parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--inVCF", type=existing_file, required=True, dest="inVCF", metavar="VCF", help="input VCF")
    parser.add_argument("-o", "--outVCF", type=existing_file, required=True, dest="outVCF", metavar="VCF", help="output VCF")
    parser.add_argument("-t", "--maxthreads", type=str, dest="maxthreads", default=DEF_MAX_THREADS, help="maximum threads [default: %(default)s]")
    parser.add_argument("-c", "--configF", type=existing_file, required=True, dest="configF", help="config file (optional)")
    parser.add_argument("-l", "--luaF", type=existing_file, required=False, dest="luaF", help="LUA file (optional)")
    args = parser.parse_args() 
    
    # load config
    config = json.load(open(args.config))
    
    # ensure dirs
    outdir = os.path.abspath(os.path.join(args.vcf, os.pardir))
    if not os.path.exists(outdir):
            print("Creating dir " + outdir)
            os.makedirs(outdir)
    log = os.path.join(outdir, "vcfanno-pipeline.log")
    success = True    

    # start 
    logging.basicConfig(filename=log, level=logging.DEBUG)                
    
    logging.info("==========================================================")
    logging.info("VCFanno pipeline")
    logging.info("==========================================================")
    logging.info("Started script at %s.", str(datetime.date.today()))
        
    runVcfanno(inVCF=args.inVCF, outVCF=args.outVCF, configF=args.configF, luaF=args.luaF, threads=args.maxthreads)
    
    # check success
    if success:     
        print ("Finished.")
        logging.info("Finished script successfully at %s.", str(datetime.date.today()))    
    else :
        print ("Pipeline had ERRORS, see log file")
        sys.exit(1)            

