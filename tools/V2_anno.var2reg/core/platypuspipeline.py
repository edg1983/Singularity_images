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

  platypus pipeline

  Copyright 2018 Niko Popitsch. All rights reserved.
  
  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0
  
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''
DEF_MAX_THREADS = 4


def runPLAT(bam, ref, vcf, roi=None, threads=DEF_MAX_THREADS, additionalParameters=[]):
    
    checkTools({"platypus", "ngs-tools", "genomic-tools", "vcf-sort" })
    
    success = True
    
    if files_exist(vcf):
        print("VCF file " + vcf + " already exists! SKIPPING re-creation...")
        logging.warn("VCF file " + vcf + " already exists! SKIPPING re-creation...")
        return success
    
    todel = []
    
    if roi is not None and roi.endswith(".gz"):
        # unzip zipped bed file
        tmp = roi + ".unziped.bed"
        success = success and pipelineStep(roi, tmp, [ "gunzip", "-c", roi], shell=True, stdout=tmp)
        roi = tmp 
        todel += [roi]
            
    # run platypus
    prefix, ext = os.path.splitext(os.path.abspath(vcf))
    iscompressed = ext.endswith('.gz')
    vcfraw = prefix + ".raw.vcf"
    vcfplatypuslog = prefix + ".plat.log"
    input = bam
    if (type(bam) is list) :
        input = ",".join(bam)
    cmd = [ getTool("platypus"),
        "--bamFiles=" + input,
        "--minFlank=0",
        "--nCPU=" + str(threads),
        "-o", vcfraw,
        "--refFile=" + ref,
        "--logFileName=" + vcfplatypuslog]
    if (roi is not None):
        cmd += ["--regions=" + roi]
    success = success and pipelineStep(None, vcfraw, cmd, shell=True)
            
    todel += [vcfraw]
    
    vcfsorted = prefix + ".sorted.vcf"    
    success = success and pipelineStep(vcfraw, vcfsorted, [ getTool("vcf-sort"), vcfraw], shell=True, stdout=vcfsorted)
    
    
    # fix CONTIG headers
    vcfsortedc = vcfsorted+".CONTIG.vcf"
    success = success and pipelineStep(vcfsorted, vcfsortedc, [ getTool("genomic-tools"),
                                                  "VCFTools", "addContigHeaders",
                                                  "-i", vcfsorted,
                                                  "-o", vcfsortedc,
                                                  "-r", ref], shell=True)
    todel += [vcfsorted]

    
    if iscompressed:
        bgzip(vcfsortedc, outFile=vcf, index=True, delinFile=True)
        todel += [vcfsortedc]
    else:
        os.rename(vcfsortedc, vcf)
    
    # create stats
    stats = prefix + ".stats"
    success = success and pipelineStep(vcf, stats, [ getTool("ngs-tools"),
                                                  "VCFTools", "stats",
                                                  "-i", vcf,
                                                  "-o", stats], shell=True)
    
    removeFile(todel)

    return success


if __name__ == "__main__":
    # ============================================================================
    parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--bam", type=existing_file, required=True, dest="bam", metavar="bam", help="input BAM")
    parser.add_argument("-o", "--vcf", type=str, required=True, dest="vcf", metavar="vcf", help="output VCF file (may be bgziped)")
    parser.add_argument("-t", "--maxthreads", type=str, dest="maxthreads", default=DEF_MAX_THREADS, help="maximum threads [default: %(default)s]")
    parser.add_argument("-r", "--ref", type=existing_file, required=True, dest="ref", help="Reference fasta")
    parser.add_argument("-c", "--config", type=existing_file, required=True, dest="config", help="Configuration file [default: %(default)s]")
    parser.add_argument("-targets", "--ROI", dest="roi", required=False, help="Regions to be considered (e.g., targeted by the used WES kit) [default: %(default)s]")
    args = parser.parse_args() 
    
    # load config
    config = json.load(open(args.config))
    
    # ensure dirs
    outdir = os.path.abspath(os.path.join(args.vcf, os.pardir))
    if not os.path.exists(outdir):
            print("Creating dir " + outdir)
            os.makedirs(outdir)
    log = os.path.join(outdir, "Platypus-pipeline.log")
    success = True    

    # start 
    logging.basicConfig(filename=log, level=logging.DEBUG)                
    
    logging.info("==========================================================")
    logging.info("PLATYPUS pipeline")
    logging.info("==========================================================")
    logging.info("Started script at %s.", str(datetime.date.today()))
        
    runPLAT(bam=args.bam, ref=args.ref, vcf=args.vcf, roi=args.roi, threads=args.maxthreads, additionalParameters=[])
    
    # check success
    if success:     
        print ("Finished.")
        logging.info("Finished script successfully at %s.", str(datetime.date.today()))    
    else :
        print ("Pipeline had ERRORS, see log file")
        sys.exit(1)            

