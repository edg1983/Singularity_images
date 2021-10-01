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

# Necessary for including python modules from a parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.utils import *
from core.config import *

#============================================================================

pipename = "SnpEff annotation pipeline"

#============================================================================
usage = pipename + '''                           

  bcltools/genomic-tools/snpEff annotaion pipeline. Takes about 90min/genome with default settings.

  Copyright 2018 CCRI. All rights reserved.
  
  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0
  
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''

REF = "/data/ref/"

#
# Additional VCF info headers
#
infos = []
infos += [vcfpy.OrderedDict([("ID", "CADD_PHRED"), ("Number", "1"), ("Type", "Float"), ("Description", "CADD phred score")])]
infos += [vcfpy.OrderedDict([("ID", "CADD_RawScore"), ("Number", "1"), ("Type", "Float"), ("Description", "CADD raw score")])]

DEF_MAX_THREADS = 16

start_time = time.time()
success = True

    
def checksuccess(stage):
    global success
    global start_time
    elapsed_time = time.time() - start_time
    if not success:
        logging.error("Pipeline failed at stage: " + stage)
        sys.exit("Pipeline failed at stage: " + stage)
    else:
        logging.info("-----------------------------------------------------------------------------")
        logging.info("Finished stage " + stage + " in " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
        logging.info("-----------------------------------------------------------------------------")
        

def runSnpEffAnno(inVcf, outVcf, genome, source, snpEffConfig=None, threads=DEF_MAX_THREADS, includeNonPass=True, overwrite=False):
    config = checkTools( {"snpEff", "snpSift", "genomic-tools", "bcftools"})
    
    if genome not in config["ref"]:
        logging.error("Genome " + genome + " not supported yet!")
        return False
    
    start_time = time.time()
    success = True
    tmp = outVcf + ".tmp.vcf"
    tmp2 = outVcf + ".tmp2.vcf"
    if overwrite or not files_exist([outVcf]):
        # ###############################################################################
        #             NORMALIZE VCF (e.g., remove non-canonical chroms, non pass variants, etc.
        # ###############################################################################
        f = []
        for c in config["ref"][genome]["chr"]:    
            f += ["( CHROM = '" + c + "' )"]
        filter = "|".join(f)
        if not includeNonPass:
            filter = "(" + filter + ") & (( na FILTER ) | (FILTER = 'PASS'))"
        
        cmd = [getTool("snpSift"), "filter", "\"" + filter + "\"", "-f", inVcf ]
        success = success and pipelineStep(inVcf, tmp, cmd, shell=True, stdout=tmp)
        checksuccess("Normalization")

        # ###############################################################################
        #             add GT info (for strelka files!)
        # ###############################################################################
        if source.lower() == "strelka":
            logging.info("Adding GT info to strelka VCF")
            cmd = [getTool("genomic-tools"), "vcftools", "strelkaAddGT", "-i", tmp, "-o", tmp2 ]
            success = success and pipelineStep(tmp, tmp2, cmd, shell=True)
            checksuccess("AddGenotypeInfo")       
            os.rename(tmp2, tmp)
        
        # ###############################################################################
        #             add cons data
        # ###############################################################################
        if (genome in config["cons"]):
            cmd = [getTool("genomic-tools"), "vcftools", "consAddPara", "-i", tmp, "-c", config["cons"][genome], "-o", tmp2, "-t", str(threads)]
            success = success and pipelineStep(tmp, tmp2, cmd, shell=True)
            checksuccess("AddConservationScores")       
            os.rename(tmp2, tmp)
        else:
            logging.warn("Skipped CONS annotation as not defined for genome " + genome)
        
        # ###############################################################################
        #             add ROI info 
        # ###############################################################################
        if (genome in config["roi"]):
            for r, rf in config["roi"][genome].items():
                cmd = [getTool("bcftools"), "annotate", "-a", rf, "-c", "CHROM,FROM,TO", "--mark-sites", "\"" + r + "\"", tmp ]
                success = success and pipelineStep(tmp, tmp2, cmd, shell=True, stdout=tmp2)
                checksuccess("Annotate " + r) 
                os.rename(tmp2, tmp)
        else:
            logging.warn("Skipped ROI annotation as not defined for genome " + genome)

        # ###############################################################################
        #             add anno info 
        # ###############################################################################
        if (genome in config["anno"]):
            for a, af in config["anno"][genome].items():
                cmd = [getTool("snpSift"), "Annotate", "-a", "-noId", "-tabix", "-a", "-name", a + "_", "-info", af[1], af[0], tmp ]
                success = success and pipelineStep(tmp, tmp2, cmd, shell=True, stdout=tmp2)
                checksuccess("Annotate " + a) 
                os.rename(tmp2, tmp)
        else:
            logging.warn("Skipped additions annotations as not defined for genome " + genome)

        # ###############################################################################
        #             add ID info 
        # ###############################################################################
        if (genome in config["known"]):
            for a, af in config["known"][genome].items():
                cmd = [getTool("snpSift"), "Annotate", "-a", "-id", "-noInfo", "-tabix", "-a", af[0], tmp ]
                success = success and pipelineStep(tmp, tmp2, cmd, shell=True, stdout=tmp2)
                checksuccess("Annotate " + a) 
                os.rename(tmp2, tmp)
        else:
            logging.warn("Skipped known variant annotations as not defined for genome " + genome)
                    
        # ###############################################################################
        #             add missing headers
        # ###############################################################################
        with vcfpy.Reader.from_path(tmp) as r:
            for i in infos:
                r.header.add_info_line(i)
            with vcfpy.Writer.from_path(tmp2, r.header) as w:
                for record in r:
                    w.write_record(record)
        os.rename(tmp2, tmp)    
        checksuccess("Add headers") 
        
        # ###############################################################################
        #             run snpEff
        # ###############################################################################
        if (genome in config["ref"]):
            cmd = [getTool("snpEff"), "ann", config["ref"][genome]["snpEff"]]
            cmd += ["-i", "vcf"]
            cmd += ["-t"]
            cmd += ["-csvStats", outVcf + ".csvStats.txt"]
            if (snpEffConfig is not None):
                cmd += ["-config", snpEffConfig]
            cmd += ["-noLog"]
            cmd += [tmp]
            success = success and pipelineStep(tmp, tmp2, cmd, shell=True, stdout=tmp2)
            checksuccess("snpEff") 
            os.rename(tmp2, tmp)
        else:
            logging.warn("Skipped snpEff annotations as not defined for genome " + genome)
        
        # ###############################################################################
        #             compress results
        # ###############################################################################
        if outVcf.endswith(".gz"):
            bgzip(tmp, outVcf, index=True, override=True, delinFile=True)
        else:
            os.rename(tmp, outVcf)
    else:
        logging.info("Skipped creation of " + outVcf + " as file already exists")

    # ###############################################################################
    #                Cleanup
    # ###############################################################################
    if (success):  
        removeFile(tmp2)
    else:
        print("Pipeline had errors!")         
    return success     
        
        
##############################################################################################
#        Commandline
##############################################################################################
if __name__ == "__main__":   
    parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--inVcf", required=True, dest="inVcf", help="path to the input VCF file [default: %(default)s]", type=existing_file, metavar="<VCF>")
    parser.add_argument("-o", "--outVcf", required=True, dest="outVcf", help="path to the output VCF file [default: %(default)s]", type=str, metavar="<VCF>")
    parser.add_argument("-g", "--genome", required=True, type=str, dest="genome", help="The used reference genome (must be in config file)")
    parser.add_argument("-s", "--source", type=str, default="platypus", dest="source", help="The used variant caller. If 'strelka', a GT field calculated from NT and SGT fields will be added to the VCF.")
    parser.add_argument("-c", "--config", type=existing_file, default=None, required=True, dest="config", help="Config file.")
    parser.add_argument("--snpEffConfig", type=existing_file, default=None, required=False, dest="snpEffConfig", help="Optional snpEff config file.")
    parser.add_argument("-l", "--log", required=False, type=str, dest="logfilename", help="logfile")
    
    parser.add_argument("--includeNonPass", required=False, action='store_true', dest="includeNonPass", help="If set, non-pass variants will be kept in the file (and annotated)")
    parser.add_argument("-t", "--maxthreads", type=str, dest="maxthreads", default=DEF_MAX_THREADS, help="maximum threads [default: %(default)s]")
    
    args = parser.parse_args()
     
    # load config and merge with default tool config to enable system-specific program locations
    global config
    config = loadUserConfig(args.config)
    
    outdir = os.path.abspath(os.path.join(args.outVcf, os.pardir))
    #============================================================================
    # ensure dirs
    #============================================================================
    if not outdir.endswith("/"):
        outdir += "/"
    if not os.path.exists(outdir):
            print("Creating dir " + outdir)
            os.makedirs(outdir)
    log = os.path.join(outdir, pipename.replace(" ", "-") + ".log")
    success = True    
    
    # start 
    logging.basicConfig(filename=log, level=logging.DEBUG)                
    logging.info("==========================================================")
    logging.info(pipename)
    logging.info("==========================================================")
    logging.info("Effective configuration:")
    logging.info(getConfig())
    logging.info("Started script at %s.", str(datetime.date.today()))
    
    logging.info("NOTE: Using " + args.genome + "annotation files!")
    
    success = success and runSnpEffAnno(inVcf=args.inVcf,
                                        outVcf=args.outVcf,
                                        genome=args.genome,
                                        source=args.source,
                                        maxthreads=args.maxthreads,
                                        includeNonPass=args.includeNonPass,
                                        snpEffConfig=args.snpEffConfig)
    
    if not success:
        sys.exit("Pipeline failed")
    else:
        elapsed_time = time.time() - start_time
        logging.info("-----------------------------------------------------------------------------")
        logging.info("Finished pipeline in " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
        logging.info("-----------------------------------------------------------------------------")
