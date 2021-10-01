'''
Created on February, 2018

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
from core.config import *


usage = '''python pipeline.py                             

  bam2fastq pipeline

  Copyright 2018 Niko Popitsch. All rights reserved.
  
  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0
  
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.
  
  Maps RNAseq reads.

USAGE
'''
DEF_MAX_THREADS = 4


def runBAM2FASTQ(bam, fq1, fq2=None, threads=DEF_MAX_THREADS):
    checkTools({"samtools", "bedtools"})
    success = True
       
    if not files_exist([fq1]):
        # convert to fastq
        outdir = os.path.abspath(os.path.join(fq1, os.pardir)) + "/"
        sns = outdir + replaceExtension(getFilename(bam), ".unmapped.namesort.bam")     
        cmd = [getTool("samtools"), "sort", "-@" + str(threads), "-n", bam, "-o", sns]
        if not files_exist([sns]):
            success = success and pipelineStep([bam], [sns], cmd, shell=True)
        
        dozip = fq1.endswith(".gz")
        if dozip:
            fq1 = replaceExtension(fq1, "")
            if fq2 is not None:
                fq2 = replaceExtension(fq2, "")
            
        cmd = [getTool("bedtools"), "bamtofastq", "-i", sns, "-fq", fq1]
        if fq2 is not None:
            cmd += ["-fq2", fq2]
        success = success and pipelineStep([sns], [fq1, fq2], cmd, shell=True)
        if dozip:
            bgzip(fq1, delinFile=True)
            if fq2 is not None:
                bgzip(fq2, delinFile=True)
        removeFile(sns)
    else:
        logging.info("Skipped bam2fastq as " + fq1 + " already exists")    
    return success


if __name__ == "__main__":
    # ============================================================================
    parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--bam", type=existing_file, required=True, dest="bam", metavar="bam", help="input BAM (unaligned)")
    parser.add_argument("-f1", "--fastq1", type=str, required=True, dest="fq1", metavar="fq1", help="1st read pair")
    parser.add_argument("-f2", "--fastq2", type=str, required=False, dest="fq2", metavar="fq2", help="2nd read pair")
    parser.add_argument("-c", "--config", type=existing_file, required=True, dest="config", help="Configuration file [default: %(default)s]")
    parser.add_argument("-t", "--maxthreads", type=str, dest="maxthreads", default=DEF_MAX_THREADS, help="maximum threads [default: %(default)s]")
    args = parser.parse_args() 
    
    # load config
    loadUserConfig(args.config)

    # ensure dirs
    outdir = os.path.abspath(os.path.join(args.fq1, os.pardir))
    if not os.path.exists(outdir):
            print("Creating dir " + outdir)
            os.makedirs(outdir)
    log = os.path.join(outdir, "BAM2FASTQ-pipeline.log")
    success = True    
     
    # start 
    logging.basicConfig(filename=log, level=logging.DEBUG)                
    
    logging.info("==========================================================")
    logging.info("BAM2FASTQ pipeline")
    logging.info("==========================================================")
    logging.info("Started script at %s.", str(datetime.date.today()))
    
    runBAM2FASTQ(args.bam, args.fq1, args.fq2, threads=args.maxthreads)

    # check success
    if success:     
        print("Finished.")
        logging.info("Finished script successfully at %s.", str(datetime.date.today()))    
    else :
        print("Pipeline had ERRORS, see log file")
        sys.exit(1)            

