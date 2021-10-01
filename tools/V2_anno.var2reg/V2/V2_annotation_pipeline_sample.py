#!/usr/bin/env python3
'''
@author: niko
'''

from argparse import ArgumentParser, RawDescriptionHelpFormatter, \
    ArgumentTypeError
import csv, datetime, time, logging, sys, os, json
from subprocess import *
import vcfpy
import re

# Necessary for including python modules from a parent directory
sys.path.append("/opt/V2_core")

from core.vcfannoPipeline import runVcfanno
from core.snpEffAnnoPipeline import runSnpEffAnno
from core.snpEffFilterPipeline import runSnpEffFilter
from core.utils import *
from core.config import *
from V2_postfilter_RD import postfilter

#============================================================================

pipename = "V2 sample annotation pipeline"

DEF_MAX_THREADS = 16
start_time = time.time()
success = True

#============================================================================
usage = pipename + '''                           

  Copyright 2019 Niko Popitsch. All rights reserved.
  
  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0
  
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''


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
    start_time = time.time()


def preprocessPlatypus(inVcf, ref, outVcf, tmpdir=None):
    global success
    logging.info("Adding CONTIG info to platypus VCF")
    vcfsorted = outVcf + ".sorted.vcf"    
    cmd = [getTool("vcf-sort")]
    if tmpdir:
        cmd += ["--temporary-directory", tmpdir]
    cmd+=[inVcf]
    pipelineStep(inVcf, vcfsorted, cmd, shell=True, stdout=vcfsorted)
    # fix CONTIG headers
    vcfsortedc = vcfsorted + ".CONTIG.vcf"
    success = success and pipelineStep(vcfsorted, None, [ getTool("genomic-tools"),
                                                  "VCFTools", "addContigHeaders",
                                                  "-i", vcfsorted,
                                                  "-o", vcfsortedc,
                                                  "-r", ref], shell=True)
    
    if not files_exist(vcfsortedc):
        # addContigHeaders will not create a new file if contig headers were already found!
        os.rename(vcfsorted, outVcf)
    else:
        os.rename(vcfsortedc, outVcf)
        removeFile(vcfsorted)  
    return outVcf

#
# preprocess deepvariant VCFs: remove ID entries
#
def preprocessDeepvariant(inVcf):
    global success
    logging.info("Removing deepvariant entries from ID field")
    vcf_no_id = inVcf + ".no_id.vcf"
    success = success and pipelineStep(inVcf, vcf_no_id, [ getTool("genomic-tools"),
                                                  "VCFTools", "moveId2InfoField",
                                                  "-n", "old_deepvariant_ids",
                                                  "-i", inVcf,
                                                  "-o", vcf_no_id], shell=True)
    
    if success:
        os.rename(vcf_no_id, inVcf)
    return inVcf

#
# preprocess deepvariant VCFs: remove ID entries
#
def preprocessGatk(inVcf):
    global success
    logging.info("Removing gatk entries from ID field")
    vcf_no_id = inVcf + ".no_id.vcf"
    success = success and pipelineStep(inVcf, vcf_no_id, [ getTool("genomic-tools"),
                                                  "VCFTools", "moveId2InfoField",
                                                  "-n", "old_gatk_ids",
                                                  "-i", inVcf,
                                                  "-o", vcf_no_id], shell=True)
    
    if success:
        os.rename(vcf_no_id, inVcf)
    return inVcf

#
# preprocess strelka VCFs: add GT field
#
def preprocessStrelka(inVcf, outVcf):
    global success
    logging.info("Adding GT info to strelka VCF")
    cmd = [getTool("genomic-tools"), "vcftools", "strelkaAddGT", "-i", inVcf, "-o", outVcf ]
    success = success and pipelineStep(tmp, tmp2, cmd, shell=True)    
    return outVcf

    
#============================================================================
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-i", "--in", type=existing_file, required=True, dest="infile", metavar="infile", help="Input sample sheet")
parser.add_argument("-o", "--outdir", type=str, required=True, dest="outdir", metavar="outdir", help="output directory")
parser.add_argument("-t", "--maxthreads", type=str, dest="maxthreads", default=DEF_MAX_THREADS, help="maximum threads [default: %(default)s]")
parser.add_argument("-c", "--config", type=existing_file, required=True, dest="config", help="Configuration file [default: %(default)s]")
parser.add_argument("-s", "--snpEffConfig", type=existing_file, required=False, dest="snpEffConfig", help="Optional snpEff configuration file [default: %(default)s]")
parser.add_argument("--includeNonPass", required=False, action="store_true", default=False, dest="includeNonPass", help="If set, non-PASS variants will not be filtered")
parser.add_argument("--noQC", required=False, action="store_true", default=False, dest="noQC", help="If set, no QC PDF will be created")
parser.add_argument("--overwrite", required=False, action="store_true", default=False, dest="overwrite", help="If set, existing files will be overwritten")

args = parser.parse_args() 
#============================================================================

# load config and merge with default tool config to enable system-specific program locations
global config
config = loadUserConfig(args.config)

# confgured linux_temp_dir?
linux_temp_dir = config["linux_temp_dir"] if "linux_temp_dir" in config else None

# ensure dirs
if not os.path.exists(args.outdir):
        print("Creating dir " + args.outdir)
        os.makedirs(args.outdir)
log = os.path.join(args.outdir, pipename.replace(" ", "-") + ".log")
if not args.outdir.endswith("/"):
    args.outdir += "/"

# start 
logging.basicConfig(filename=log, level=logging.DEBUG)                
logging.info("==========================================================")
logging.info(pipename)
logging.info("==========================================================")
logging.info("Effective configuration:")
logging.info(getConfig())
logging.info("Started script at %s.", str(datetime.date.today()))

with open(args.infile, 'rt') as tsvin:
    tsvin = csv.reader(tsvin, delimiter='\t')
    for row in tsvin:
        if (row[0].startswith("#")):
            continue
        PID = row[0]
        genome = row[1]
        source = row[2]
        inVcf = row[3]      
        outdir = args.outdir + PID 
        if not files_exist(inVcf):
            print("Skipping configured input vcf as file not found: " + inVcf)
            continue
        if not files_exist(outdir):
            os.makedirs(outdir)
        tmpdir = outdir + "/tmp"
        if not files_exist(tmpdir):
            os.makedirs(tmpdir)
        qcdir = outdir + "/qc"
        if not files_exist(qcdir):
            os.makedirs(qcdir)
        ds = PID + "." + genome

        vcfa = tmpdir + "/" + ds + ".anno.vcf.gz"  # temporary unfiltered VCF
        vcfaf = outdir + "/" + ds + ".anno+fil.vcf.gz"
        tsv = tmpdir + "/" + ds + ".anno+fil.tsv"  # temporary snpSift output
        final = outdir + "/" + ds + ".final.tsv"
        qcpdf = qcdir + "/" + ds + ".qc.pdf"
        
        logging.info("-------------------------")
        logging.info(ds)
        logging.info("-------------------------")
        print(ds)

        if args.overwrite or not files_exist([vcfa]):
            tmp = tmpdir + "/" + ds + ".tmp1.vcf"
            tmp2 = tmpdir + "/" + ds + ".tmp2.vcf"
            
            # ###############################################################################
            #             NORMALIZE VCF (e.g., remove non-canonical chroms, non pass variants, etc.)
            # ###############################################################################
            f = []
            for c in config["ref"][genome]["chr"]:    
                f += ["( CHROM = '" + c + "' )"]
            filter = "|".join(f)
            if not args.includeNonPass:
                filter = "(" + filter + ") & (( na FILTER ) | (FILTER = 'PASS'))"
            
            cmd = [getTool("snpSift"), "filter", "\"" + filter + "\"", "-f", inVcf ]
            success = success and pipelineStep(inVcf, tmp, cmd, shell=True, stdout=tmp)
            checksuccess("Normalization")
    
            # ###############################################################################
            #             PREPROCESS VCF (e.g., add GT info to somatic VCF)
            # ###############################################################################
            if source.lower() == "platypus":
                # ###############################################################################
                # fix CONTIG headers               
                # ###############################################################################
                preprocessPlatypus(inVcf=tmp, ref=config["ref"][genome]["FASTA"], outVcf=tmp, tmpdir=linux_temp_dir)
            elif source.lower() == "deepvariant":
                # ###############################################################################
                # remove IDs as they are not informative        
                # ###############################################################################
                preprocessDeepvariant(inVcf=tmp)
            elif source.lower() == "gatk":
                # ###############################################################################
                # remove IDs as they are not informative        
                # ###############################################################################
                preprocessGatk(inVcf=tmp)
            elif source.lower() == "strelka":
                # ###############################################################################
                #                add GT info (for strelka files!)
                # ###############################################################################
                preprocessStrelka(inVcf=tmp, outVcf=tmp)        
            checksuccess("Preprocessing")
            
            # ###############################################################################
            #             Annotate with snpEFF
            # ###############################################################################
            if (genome in config["ref"]):
                cmd = [getTool("snpEff"), "ann", config["ref"][genome]["snpEff"]]
                cmd += ["-i", "vcf"]
                cmd += ["-t"]
                cmd += ["-csvStats", tmpdir + "/" + ds + ".snpEff.ann.csvStats.txt"]
                if (args.snpEffConfig is not None):
                    cmd += ["-config", args.snpEffConfig]
                cmd += ["-noLog"]
                cmd += [tmp]
                success = success and pipelineStep(tmp, tmp2, cmd, shell=True, stdout=tmp2)
                checksuccess("snpEff") 
                os.rename(tmp2, tmp)
            else:
                logging.warn("Skipped snpEff annotations as not defined for genome " + genome)
    
            # ###############################################################################
            #             Annotate with vcfanno
            # ###############################################################################  
            vcfaconf = tmpdir + "/" + ds + ".vcfanno.config"
            vcfalua = tmpdir + "/" + ds + ".vcfanno.lua"
            # create lua file
            with open(vcfalua, "w") as f: 
                f.write("function setid(pre,...)\n")
                f.write(" local t = {...}\n")
                f.write(" local res = {}\n")
                f.write(" local seen = {}\n")
                f.write(" for idx, ids in pairs(t) do\n")
                f.write("  local sep=\",\"\n")  # split comma-separated IDs
                f.write("  for v in string.gmatch(ids, \"([^\"..sep..\"]+)\") do\n")
                f.write("   for i, v in pairs(t) do\n")
                f.write("    if v ~= \".\" and v ~= nil and v ~= \"\" then\n")
                f.write("     if seen[v] == nil then\n")
                f.write("      res[#res+1] = string.gsub(pre[i] .. v, \",\", \";\")\n")
                f.write("      seen[v] = true\n")
                f.write("     end\n")
                f.write("    end\n")
                f.write("   end\n")
                f.write("  end\n")
                f.write(" end\n")
                f.write(" return table.concat(res, \";\")\n")
                f.write("end\n")
            logging.info("Created LUA file at " + vcfalua)    
            
            # create vcfanno config file     
            with open(vcfaconf, "w") as f: 
                if ("roi" in config and genome in config["roi"]):
                    for r, rf in config["roi"][genome].items():
                        f.write("[[annotation]]\n")
                        f.write("file=\"" + rf + "\"\n")
                        f.write("columns = [3]\n")  # use column 3 as some BED do not contain 4 columns
                        f.write("ops = [\"flag\"]\n")
                        f.write("names = [\"" + r + "\"]\n")
                        f.write("\n")
                else:
                    logging.warn("Skipped ROI annotation as not defined for genome " + genome)
                    
                if ("af" in config and genome in config["af"]):
                    file2block = {}
                    postanno = {}
                    for id in config["af"][genome].keys():
                        for a, af in config["af"][genome][id].items():
                                                       
                            block = file2block[af[0]] if af[0] in file2block else [[], [], []]  # fields, ops, names
                            
                            if af[1].count(",") == 1:  # configured AC,AN. We need a postannoblock
                                postanno[id + "_" + a]=[]
                                for x in af[1].split(","):
                                    postanno[id + "_" + a].append( "\"" + id + "_" + a + "_" + x +"\"" )
                                
                            for x in af[1].split(","):
                                block[0].append("\"{0}\"".format(x))
                            if (len(af) > 2):
                                for x in af[1].split(","):
                                    block[1].append("\"{0}\"".format(x))
                            else:
                                for x in af[1].split(","):
                                    block[1].append("\"max\"")
                            for x in af[1].split(","):
                                block[2].append("\"" + id + "_" + a + "_{0}\"".format(x))
                            file2block[af[0]] = block
                    # write blocks
                    for file in file2block:    
                        block = file2block[file]
                        f.write("[[annotation]]\n")
                        f.write("file=\"" + file + "\"\n")
                        # select AF INFO fields
                        f.write("fields = [" + ", ".join(block[0]) + "]\n")
                        f.write("ops = [" + ", ".join(block[1]) + "]\n")
                        f.write("names = [" + ", ".join(block[2]) + "]\n")
                        f.write("\n")
                        
                    # write postanno blocks
                    for id in postanno:
                        # divide ac/an
                        f.write("[[postannotation]]\n")
                        f.write("fields=[" + ", ".join(postanno[id]) + "]\n")
                        f.write("name=\"" + id + "_AF\"\n")
                        f.write("op=\"div2\"\n")
                        f.write("type=\"Float\"\n")
                        f.write("\n")
                else:
                    logging.warn("Skipped AF annotation as not defined for genome " + genome)
                    
                if ("anno" in config and genome in config["anno"]):
                    for a, af in config["anno"][genome].items():
                        if af[0].endswith("bed") or af[0].endswith("bed.gz"):
                            # BED
                            f.write("[[annotation]]\n")
                            f.write("file=\"" + af[0] + "\"\n")
                            f.write("columns = [4]\n")  # copy value from 4th column (=BED name)
                            f.write("ops = [\"self\"]\n")
                            f.write("names = [\"" + a + "\"]\n")
                            f.write("\n")
                        elif af[0].endswith("tsv") or af[0].endswith("tsv.gz") or af[0].endswith("txt") or af[0].endswith("txt.gz"):
                            # TSV: 0=file, 1=columns, 2=names, 3=ops
                            f.write("[[annotation]]\n")
                            f.write("file=\"" + af[0] + "\"\n")
                            # columns indices
                            f.write("columns = [" + ", ".join("{0}".format(e) for e in af[1].split(",")) + "]\n")
                            
                            # column names
                            if (len(af) > 2):
                                f.write("names = [" + ", ".join("\"" + a + "_{0}\"".format(e) for e in af[2].split(",")) + "]\n")
                            else: 
                                # use column indices as names...
                                f.write("names = [" + ", ".join("\"" + a + "_{0}\"".format(e) for e in af[1].split(",")) + "]\n")
                            # ops
                            if (len(af) > 3):
                                f.write("ops = [" + ", ".join("\"{0}\"".format(e) for e in af[3].split(",")) + "]\n")
                            else:
                                f.write("ops = [" + ", ".join("\"self\"".format(e) for e in af[1].split(",")) + "]\n")
                            f.write("\n")
                        else:
                            # VCF
                            f.write("[[annotation]]\n")
                            f.write("file=\"" + af[0] + "\"\n")
                           # select INFO fields
                            f.write("fields = [" + ", ".join("\"{0}\"".format(e) for e in af[1].split(",")) + "]\n")
                            if (len(af) > 2):
                                f.write("ops = [" + ", ".join("\"{0}\"".format(e) for e in af[2].split(",")) + "]\n")
                            else:
                                f.write("ops = [" + ", ".join("\"self\"".format(e) for e in af[1].split(",")) + "]\n")
                            f.write("names = [" + ", ".join("\"" + a + "_{0}\"".format(e) for e in af[1].split(",")) + "]\n")
                            f.write("\n")
                else:
                    logging.warn("Skipped ANNO annotation as not defined for genome " + genome)
            
                if ("known" in config and genome in config["known"]):
                    knownFields=[]
                    knownFieldsQuote=[]
                    prefixes=[]
                    for a, af in config["known"][genome].items():
                        # annotate =INFO field
                        f.write("[[annotation]]\n")
                        f.write("file=\"" + af[0] + "\"\n")
                        f.write("fields=[\"ID\"]\n")
                        f.write("names=[\"" + a + "_ID\"]\n")
                        f.write("ops=[\"self\"]\n")
                        knownFieldsQuote+=["\"" + a + "_ID\""]
                        knownFields+=[a + "_ID"]
                        prefix = ""
                        if (len(af) > 1):
                            prefix = af[1]
                        prefixes+=["'" + prefix +"'"]
            
                    # postannotate: move from INFO fields to ID
                    f.write("[[postannotation]]\n")
                    f.write("name=\"ID\"\n")
                    f.write("fields=[" + ",".join(knownFieldsQuote)+"]\n")
                    #f.write("op=\"lua:setid(" + prefixes + "', " + a + "_ID)\"\n")
                    f.write("op=\"lua:setid({"+",".join(prefixes)+"}, "+",".join(knownFields)+")\"\n")
                    f.write("type=\"String\"\n")

                    # postannotate: delete INFO fields 
                    f.write("[[postannotation]]\n")
                    f.write("fields=[" + ",".join(knownFieldsQuote)+"]\n")
                    f.write("op=\"delete\"\n")
                    f.write("\n")
            logging.info("Created vcfanno configuration file at " + vcfaconf)   
             
            # run vcfanno
            runVcfanno(inVCF=tmp, outVCF=tmp2, configF=vcfaconf, luaF=vcfalua, threads=args.maxthreads, settings=["GOGC=1000", "IRELATE_MAX_CHUNK=8000", "IRELATE_MAX_GAP=1000"])
            checksuccess("anno") 
            os.rename(tmp2, tmp)
            
    #         # ###############################################################################
    #         #             add missing headers
    #         # ###############################################################################
    #         with vcfpy.Reader.from_path(tmp) as r:
    #             for i in infos:
    #                 r.header.add_info_line(i)
    #             with vcfpy.Writer.from_path(tmp2, r.header) as w:
    #                 for record in r:
    #                     w.write_record(record)
    #         os.rename(tmp2, tmp)    
    #         checksuccess("Add headers") 
    
            # ###############################################################################
            #             compress results
            # ###############################################################################
            if vcfa.endswith(".gz"):
                bgzip(tmp, vcfa, index=True, override=True, delinFile=True, maxthreads=args.maxthreads)
            else:
                os.rename(tmp, vcfa)
        else:
            print("Annotated VCF exists, will not re-create: " + vcfa)

        if args.overwrite or not files_exist([vcfaf]):
            # ###############################################################################
            #             Filter VCF (e.g., remove non-canonical chroms, non pass variants, etc.
            # ###############################################################################
            snpEffFilterF = tmpdir + "/" + ds + ".snpEffFilter.txt"
            with open(snpEffFilterF, "w") as f: 
                f.write(config["ref"][genome]["filter"])
            tmp = vcfaf + ".tmp.vcf"
            cmd = [getTool("snpSift"), "filter"]
            cmd += ["-f", vcfa ]
            cmd += ["-e", snpEffFilterF]
            cmd += ["-i", "V2filter"]
            cmd += ["-p"]
            success = success and pipelineStep([vcfa, snpEffFilterF], tmp, cmd, shell=True, stdout=tmp)
            if vcfaf.endswith(".gz"):
                bgzip(tmp, vcfaf, index=True, override=True, delinFile=True, maxthreads=args.maxthreads)
            else:
                os.rename(tmp, vcfaf)
            checksuccess("Filtering")
    
        reader = vcfpy.Reader.from_path(vcfaf)
        samples = reader.header.samples.names
        print("Samples: %s" % samples)
                    
        if args.overwrite or not files_exist([tsv + ".gz"]):
            # ###############################################################################
            #             Extract fields to TSV
            # ###############################################################################
            # "ANN[*].RANK", "ANN[*].CDNA_POS", "ANN[*].CDNA_LEN", "ANN[*].CDS_POS", "ANN[*].CDS_LEN","ANN[*].AA_POS", "ANN[*].AA_LEN", "ANN[*].DISTANCE",
            headers = []
            for h in ["CHROM", "POS", "ID", "REF", "ALT", "FILTER", "AF", "AC", "DP", "MQ",
                      "ANN[*].ALLELE", "ANN[*].EFFECT", "ANN[*].IMPACT", "ANN[*].GENE", "ANN[*].GENEID",
                      "ANN[*].FEATURE", "ANN[*].FEATUREID", "ANN[*].BIOTYPE", "ANN[*].HGVS_C",
                      "ANN[*].HGVS_P", "ANN[*].ERRORS", "GEN[*].GT", "GEN[*].GQ"]: 
                headers += ["\"" + h + "\""]
                  
            # add INFO fields from VCF

            pattern = re.compile("^[A-Za-z_][0-9A-Za-z_.]*\Z")
            for h in reader.header.info_ids():
                # check compatibility
                if pattern.match(h):
                    headers += ["\"" + h + "\""]
                else:
                    logging.warn("Cannot extract INFO field " + h + " as not VCF compliant")
                
            cmd = [getTool("snpSift"), "extractFields"]
            cmd += ["-s", "\",\""]
            cmd += ["-e", "\".\""]
            cmd += [vcfaf]
            cmd += headers
            success = success and pipelineStep(vcfaf, tsv, cmd, shell=True, stdout=tsv)
            bgzip(tsv, delinFile=True, override=True, maxthreads=args.maxthreads)
            checksuccess("Extracting")
 
        #
        # postprocess
        #
        if args.overwrite or not files_exist([final + ".gz"]): 

            postfilter(tsv + ".gz", final, genome, args.config, samples, overwrite=True, maxthreads=args.maxthreads)
            bgzip(final, delinFile=True, override=True, maxthreads=args.maxthreads)
        
        #
        # QC
        #
        if not args.noQC and (args.overwrite or not files_exist([qcpdf])): 
                # run QC
                rscript = os.path.dirname(os.path.realpath(__file__)) + "/V2_qc.R"
                cmd = [ "Rscript", rscript, final + ".gz", args.config, genome, qcpdf]
                pipelineStep(None, None, cmd, shell=True)
                checksuccess("QC")
                        
# check success
if success:     
    print ("Finished.")
    logging.info("Finished script successfully at %s.", str(datetime.date.today()))    
else :
    print ("Pipeline had ERRORS, see log file")
    sys.exit(1)            

