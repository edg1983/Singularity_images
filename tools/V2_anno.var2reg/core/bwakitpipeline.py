'''
@author: niko
'''

from argparse import ArgumentParser, RawDescriptionHelpFormatter, \
    ArgumentTypeError
from subprocess import *
import datetime, time
import logging
import sys, os, json

# Necessary for including python modules from a parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.utils import *
from core.config import *

usage = '''python pipeline.py                             

  Human BWAKIT pipeline

  Copyright 2018 Niko Popitsch. All rights reserved.
  
  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0
  
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.
  
  Maps RNAseq reads.

USAGE
'''
DEF_MAX_THREADS = 4


def runBWAKIT(bam, ref, reads1=[], reads2=[], threads=DEF_MAX_THREADS, runFastqc=True, additionalParameters=[], id=None):
    success = True
    
    checkTools({"bwakit", "ngs-tools"})
    
    if files_exist(bam):
        print("BAM file " + bam + " already exists! SKIPPING re-creation...")
        logging.warn("BAM file " + bam + " already exists! SKIPPING re-creation...")
        return success
    
    # check input params
    if len(reads1) == 0 or len(reads1) != len(reads2):
        return False

    if not id:
        id = getFilenameNoExt(bam)
    
    RG = "\"@RG\tID:" + id + "\tSM:" + id + "\tPL:ILLUMINA\""
    
    # PAIRED END
    if (not files_exist([bam])):   
        tomerge = []
        for idx, val in enumerate(reads1):
            bamname = os.path.splitext(bam)[0] + "." + str(idx)
            bwaout = bamname + ".aln.bam"
            sorted = bamname + ".bam"
            cmd = [getTool("bwakit")]
            cmd += ["-o", bamname]
            cmd += ["-R", RG]
            cmd += ["-S"]
            cmd += ["-d"]  # deduplication with samblaster
            cmd += ["-t", str(threads)]
            cmd += [ref]
            cmd += [reads1[idx], reads2[idx]]
            cmd += ["|", "sh"]
            if (additionalParameters):
                cmd += additionalParameters
            success = success and pipelineStep([reads1, reads2], bwaout, cmd, shell=True) 
            success = success and sambamba2bam(bwaout, sorted, 
                                          sort=True,
                                          index=True,
                                          delinFile=True,
                                          maxthreads=threads) 
            tomerge += [sorted]
            
        # merge files
        merged = os.path.splitext(bam)[0] + ".merged.bam"
        if len(tomerge) == 1:
            os.rename(tomerge[0], merged)
            os.rename(tomerge[0] + ".bai", merged + ".bai")
        else:
            logging.info("Merging from multiple input files")
            success = success and sambambamerge(merged, inFiles=tomerge, maxthreads=threads) 
        
        os.rename(merged, bam)
        os.rename(merged + ".bai", bam + ".bai")
    
        flagstat(bam)
        bamstats(bam)
        
        if runFastqc:
            logging.info("Running FASTQC")
            runFASTQC(bam, os.path.dirname(bam) + "/fastqc", threads=threads)
        
    else:
        logging.info("Skipped creation of BAM file %s as it already exists.", bam) 
    return success


if __name__ == "__main__":
    # ============================================================================
    parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-r", "--ref", type=existing_file, required=True, dest="ref", help="Reference fasta")
    parser.add_argument("-i", "--id", type=str, required=False, dest="id", metavar="id", help="output id (used for filename)")
    parser.add_argument("-o", "--outdir", type=str, required=True, dest="outdir", metavar="outdir", help="output directory")
    parser.add_argument("-f1", "--fastq1", type=str, required=True, dest="fastq1", metavar="fastq1", help="1st read pair")
    parser.add_argument("-f2", "--fastq2", type=str, required=False, dest="fastq2", metavar="fastq2", help="2nd read pair")
    parser.add_argument("-t", "--maxthreads", type=str, dest="maxthreads", default=DEF_MAX_THREADS, help="maximum threads [default: %(default)s]")
    parser.add_argument("-c", "--config", type=existing_file, required=True, dest="config", help="Configuration file [default: %(default)s]")
    parser.add_argument("-nofastqc", required=False, action="store_true", default=False, dest="nofastqc", help="If set, no FASTQC analysis will be run")
    args = parser.parse_args() 
    
    
    # load config
    loadUserConfig(args.config)

    # ensure dirs
    if not os.path.exists(args.outdir):
            print("Creating dir " + args.outdir)
            os.makedirs(args.outdir)
    log = os.path.join(args.outdir, "BWAKIT-pipeline.log")
    if not args.outdir.endswith("/"):
        args.outdir += "/"
    
    success = True    
    bam = args.outdir + replaceExtension(getFilename(args.fastq1), ".bam")    
    if args.id is not None:
        bam = args.outdir + args.id + ".bam"    
      
    # start 
    logging.basicConfig(filename=log, level=logging.DEBUG)                
    
    logging.info("==========================================================")
    logging.info("BWAKIT pipeline")
    logging.info("==========================================================")
    logging.info("Started script at %s.", str(datetime.date.today()))
    
    logging.info("Mapping")
    runBWAKIT(bam=bam, ref=args.ref, reads1=[args.fastq1], reads2=[args.fastq2], threads=args.maxthreads, runFastqc=(not args.nofastqc))

    # check success
    if success:     
        print ("Finished.")
        logging.info("Finished script successfully at %s.", str(datetime.date.today()))    
    else :
        print ("Pipeline had ERRORS, see log file")
        sys.exit(1)            

