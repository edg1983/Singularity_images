#!/usr/bin/env python3
'''
@author: niko
'''

# --ids 006Bra001
# --conf /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/var2reg/var2reg.config.json
# --ped /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/noncoding/ped/{id}.ped
# --gad /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/noncoding/gado/{id}.txt
# --exo /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/noncoding/exomiser/{id}_WGS_ALL_{IM}.genes.tsv
# --reg /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/var2reg/testdata/006Bra001_TUBA1A/GRCh38_regulatory_regions.clean.006Bra001_TUBA1A.tsv.gz
# --snv /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/var2reg/testdata/006Bra001_TUBA1A/HICF2.GRCh38.anno+fil.006Bra001_TUBA1A.vcf.gz
# --cnv /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/noncoding/cnv/cnv.vcf.gz
# --gff /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/var2reg/testdata/006Bra001_TUBA1A/gencode.v32.annotation.sorted.006Bra001_TUBA1A.gff3.gz
# --gen /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/noncoding/gene_anno/gene_anno.tsv.gz
# --out /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/var2reg/results
# --include_na_models

# -ids 036EBV001
# -conf /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/var2reg/var2reg.config.json
# -ped /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/noncoding/ped/{id}.ped
# -gad /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/noncoding/gado/{id}.txt
# -exo /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/noncoding/exomiser/{id}_WGS_ALL_{IM}.genes.tsv.col12
# -reg /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/var2reg/testdata/036EBV001_DOCK8/GRCh38_regulatory_regions.clean.036EBV001_DOCK8.tsv.gz
# -snv /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/var2reg/testdata/036EBV001_DOCK8/HICF2.GRCh38.anno+fil.036EBV001_DOCK8.vcf.gz
# -cnv /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/noncoding/cnv/cnv.vcf.gz
# -gff /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/var2reg/testdata/036EBV001_DOCK8/gencode.v32.annotation.sorted.036EBV001_DOCK8.gff3.gz
# -gen /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/noncoding/gene_anno/gene_anno.tsv.gz
# -ali /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/noncoding/alias/alias2official_simple.tsv
# -o /Users/niko.popitsch/Desktop/data/projects/WTCHG/HICF2/anno/var2reg/results
# --include_na_models


from argparse import ArgumentParser, RawDescriptionHelpFormatter, \
    ArgumentTypeError
import collections
import csv, datetime, time, logging, sys, os, json
from enum import Enum
import gzip
import re
import math
import copy
from subprocess import *
import itertools
from gffutils.iterators import DataIterator
from numpy.distutils.fcompiler import none
from tqdm import tqdm
from multiprocessing import Pool, Process, Manager
from vcfpy import *
import vcfpy
import pandas as pd
import pyranges as pr
import warnings
import chunk

# Necessary for including python modules from a parent directory
sys.path.append("/opt/V2_core")
from core.utils import *

#============================================================================

pipename = "V2 var2reg pipeline"

DEF_MAX_THREADS = 16
start_time = time.time()
success = True


#============================================================================
usage = pipename + '''                           

  Copyright 2020 Niko Popitsch. All rights reserved.
  
  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0
  
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''

logo = """
 _  _  ____    _  _    _  _   __   ____  ____  ____  ____  ___
/ )( \(___ \  (_)(_)  / )( \ / _\ (  _ \(___ \(  _ \(  __)/ __) 
\ \/ / / __/   _  _   \ \/ //    \ )   / / __/ )   / ) _)( (_ \ 
 \__/ (____)  (_)(_)   \__/ \_/\_/(__\_)(____)(__\_)(____)\___/ v 0.1
=====================================================================
parallel version
"""


# @see https://stackoverflow.com/questions/2187269/print-only-the-message-on-warnings
def custom_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return str(msg) + '\n'
warnings.formatwarning = custom_formatwarning

# ==========================================================================
#    Variant filter tags
#    AF = maxPopAF too high
#    NC = no call in any affected sample
#    NI = no supported inheritance model (only if include_na_models is False)
#    LS = low score
#    ID = Incomplete data (CNVs only)
#    UT = unspecific varint type (CNVs only)
#    NP = no_pass (these variants are not included in the debug log!)
# ==========================================================================

filter_tags = ["AF", "NC", "NI", "LS", "ID", "UT", "NP" ]

# ==========================================================================
#    Scoring
# ==========================================================================
def getMax(list):
    if list is None:
        return None
    if type(list) is int:
        return list
    if type(list) is float:
        return list
    if type(list) is str:
        return None
    if list is "NA":
        return None
    if not list:
        return None
    vmax = None
    for v in list:
        if v is not None and v is not "NA":
            if vmax is None:
                try:
                    vmax = float(v)
                except ValueError:
                    vmax = None
            else:
                try:
                    vmax = float(max(float(v), vmax)) 
                except ValueError:
                    vmax = None
    return vmax

def toFloat(s):
    if s is None:
        return None
    if type(s) is float:
        return s
    if type(s) is int:
        return float(s)
    return getMax(s.split(","))

def fmt(v):
    if v is None:
        return "NA"
    return str(v)


def check_config( c ):
    # check existing sections
    for t in ["dataset_name", "input_data", "filters", "output_fields", "d_score_calc", "def_a_score", "so_term"]:
        if t not in c:
            print("Error: config section ", t, "was not found! Exiting...")
            sys.exit(1)
    # check output_fields sections
    for t in ["included_info_fields", "included_info_factors", "included_info_flags"]:
        if t not in c["output_fields"]:
            print("Error: config section ", t, "was not found! Exiting...")
            sys.exit(1)   
    # input data
    for t in ["pedigrees", "snv_vcf", "cnv_vcf", "ped_pattern", "gado_pattern","exomiser_pattern", "reg_db", "gene_gff"]:
        if t not in c["input_data"]:
            print("Error: input data section ", t, "was not found! Exiting...")
            sys.exit(1)    
    if len(c["input_data"]["pedigrees"])<=0:
            print("Error: no pedigree id configured! Exiting...")
            sys.exit(1)  
    # check for existing input files
    for t in c["input_data"]:
        if t in ["snv_vcf","cnv_vcf","alias_table","gene_anno_table","reg_db","gene_anno"]:
            checkFile(c["input_data"][t])


# calculate d_score
def calc_score(rec, conf, type, def_value):
    scores=[]
    for i,f in enumerate(conf["d_score_calc"][type]["fields"]):
        if f.startswith("score_"):
            score = calc_score(rec, conf, f[6:], def_value)
        else:
            score = rec.INFO.get(f, None) 
        if score is not None:
            scores+=[toFloat(score) / conf["d_score_calc"][type]["norm"][i]]
    if len(scores)==0:
        return def_value
    #summarize
    if conf["d_score_calc"][type]["summarisation"] == "max":
        return getMax(scores+[def_value])
    elif conf["d_score_calc"][type]["summarisation"] == "mean":
        sum=0
        for s in scores:
            sum+=s 
        return float(sum)/float(len(scores))
    sys.exit("Unknown summarisation method " + conf["d_score_calc"][type]["summarisation"])

# ==========================================================================
#    Util
# ==========================================================================
def serialize(rec):
    return rec.serialize()

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

# show progress
def show_prog(q, max_value):
    prog = tqdm(total=max_value, unit=' vars', desc='Analyzing variants')
    while 1:
        try:
            to_add = q.get(timeout=1)
            prog.n += to_add
            prog.update(0)
            q.task_done()
            if prog.n >= max_value:
                break
        except:
            continue
        
# ==========================================================================
#    SO term
# ==========================================================================
# FIXME read from JSPN config
class SoTerm():

    def __init__(self, term, d_score, a_score, type, so_id, description):
        self.term = term
        self.d_score = d_score
        self.a_score = a_score
        self.type = type
        self.so_id = so_id
        self.description = description
        
    def __repr__(self):
        return self.__str__()  
    def __str__(self):
        ret = ("%s [d:%f, a:%f]" % (self.term, self.d_score, self.a_score) )
        return (ret)   


# ==========================================================================
#    Pedigree
# ==========================================================================
class Person():

    def __init__(self, family, id, sex, affected, mum, dad):
        self.family = family
        self.id = id
        self.sex = sex
        self.affected = affected
        self.mum = mum
        self.dad = dad
        self.children = []

    def is_parent(self):
        return len(self.children) > 0

    def has_parent(self):
        return self.dad or self.mum

    def __repr__(self):
        return self.__str__()  

    def __str__(self):
        ret = ("[%s" % self.id)
        if self.mum:
            ret += " ,mum=%s" % self.mum.id
        if self.dad:
            ret += " ,dad=%s" % self.dad.id
        ret += "]"
        return (ret)   

    
class Family():
    
    def __init__(self, id):
        self.id = id
        self.members = dict()

    def affected_members(self):
        return [x for x in self.members.values() if x.affected == True]

    def unaffected_members(self):
        return [x for x in self.members.values() if x.affected == False]

    def __repr__(self):
        return self.__str__()  

    def __str__(self):
        return ('[%s: %s]' % (self.id, ", ".join(str(x) for x in self.members.values())))
    
# exomiser inheritance models
EXOMISER_IM = ["AR", "AD", "XR", "XD"]

    
class Pedigree():

    def __init__(self, pedF):
        super().__init__()
        self.name = getFilenameNoExt(pedF)
        self.families = dict()
        self.readFromFile(pedF)
        
    def readFromFile(self, pedF):
        try:
            with open(pedF, 'r') as csvfile:
                reader = csv.reader(csvfile, delimiter='\t')
                for row in reader:
                    fid = row[0]
                    pid = row[1]
                    sex = row[4]
                    aff = row[5] == '2'
                    f = self.families.get(fid, Family(fid))       
                    p = Person(f, pid, sex, aff, None, None)
                    f.members[pid] = p
                    self.families[fid] = f
        except IndexError:
            sys.exit("Error parsing pedigree " + pedF)
        # update parent/child        
        with open(pedF, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            for row in reader:
                fid = row[0]
                pid = row[1]
                did = row[2]
                mid = row[3]
                f = self.families.get(fid)       
                p = f.members.get(pid)
                if mid in f.members.keys():
                    mum = f.members.get(mid)
                    p.mum = mum
                    mum.children.append(p)
                if did in f.members.keys():
                    dad = f.members.get(did)
                    p.dad = dad
                    dad.children.append(p)

    def all_ids(self):
        ret = []
        for f in self.families.values():
            for id in f.members.keys():
                ret = ret + [id]
        return (ret) 

    def affected_ids(self):
        ret = []
        for f in self.families.values():
            for m in f.affected_members():
                ret = ret + [m.id]
        return (ret) 
    
    def get_person(self, id):
        for f in self.families.values():
            if id in f.members.keys():
                return f.members[id]
        return None

    # calculate support (=#individuals/GT that support the model-#individuals/GT that contradict the model) for all supported inheritance models.  
    # The following rules are applied:
    # * recessive: 
    #    -1    for all unaffected sample containing a HOM call
    #    +1    for all affected samples with HOM calls that are inherited from mum & dad 
    #    0     for all other samples
    # * dominant:
    #    -1    for all unaffected samples with a HOM/HET call
    #    -1    for all affected samples containing a HOM call 
    #    +1    for all affected samples with a HET call that was inherited from an affected sample if available
    #    0     for all other samples
    # * de_novo:
    #    -1    for all unaffected samples with a HET/HOM call
    #    -1    for all affected samples that inherit a call 
    #    +1    for all affected samples containing a HET call with high GQ that was not inherited from mum&dad 
    #     0    for all other samples
    # NOTE that for missing data (including calls with low GQ) we assume that GT that supports the respective inheritance model
    # independent of each other which may lead to the situation that different genotypes for the alleles are assumed per inh model.
    def calc_inh_support(self, rec, min_gq, min_gq_dnm):
        sup_dom = 0
        sup_rec = 0
        sup_dnm = 0
        sup_com = 0
        genotypes=[]
        genotype_quals=[]
        num_hq_calls = 0
        hom_aff = 0
        hom_unaff = 0
        het_aff = 0
        het_unaff = 0
        sup_comhet_mum=[]
        sup_comhet_dad=[]
        
        for pid in self.all_ids():
            p = self.get_person(pid)
            if pid not in rec.call_for_sample: # no data for this sample
                gt = None
                gq = None
            else:
                call = rec.call_for_sample[pid]
                # get genotypes
                gt = call.gt_type # hom_ref = 0, het = 1, hom_alt = 2 (which alt is untracked), uncalled = None
                gq = getMax(call.data["GQ"])
            
            # calculate VAF if possible (FIXME: for 1st allele only)
            #vaf = ad[1] / dp if ( ad is not None and len(ad)>0 and dp is not None and dp>0 ) else None
            
            genotypes.append(str(gt) if gt is not None else "NA")
            genotype_quals.append(str(gq) if gq is not None else "NA")
            
            good_gq = gq is not None and gq > min_gq
            gt_m = rec.call_for_sample[p.mum.id].gt_type if p.mum and p.mum.id in rec.call_for_sample else None
            gt_d = rec.call_for_sample[p.dad.id].gt_type if p.dad and p.dad.id in rec.call_for_sample else None
             
            is_inherited = ( gt_m is not None and gt_m > 0 ) or ( gt_d is not None and gt_d > 0 )
            is_affected = pid in self.affected_ids()

            if not good_gq or gt is None:
                continue # ignore calls with low GQ 
            

            
            # count calls
            num_hq_calls +=1
            if gt == 1:
                if is_affected:
                    het_aff += 1
                else:
                    het_unaff += 1
            if gt == 2:
                if is_affected:
                    hom_aff += 1
                else:
                    hom_unaff += 1
                    
            # calc support for comphet (only if complete genotypes!)
            if gt == 1 and is_affected:
                gq_m = getMax(rec.call_for_sample[p.mum.id].data["GQ"] if p.mum and p.mum.id in rec.call_for_sample else None)
                gq_d = getMax(rec.call_for_sample[p.dad.id].data["GQ"] if p.dad and p.dad.id in rec.call_for_sample else None)
                if gq_m is not None and gq_m > min_gq and gq_d is not None and gq_d > min_gq:
                    from_mum = False
                    from_dad = False                
                    if gt_m == 1 and gt_d == 0:
                        sup_comhet_mum.append(pid)               
                    if gt_m == 0 and gt_d == 1:
                        sup_comhet_dad.append(pid)               
                                        
            # calc supp
            if is_affected:
                if gt==0: # homref. reduces support for dom/recessive inheritance
                    sup_dnm += 0
                    sup_dom -= 1
                    sup_rec -= 1
                elif gt==1: # het
                    sup_rec -= 1
                    if gt_m is not None and gt_d is not None:
                        # check gq of parent calls
                        gq_m = getMax(rec.call_for_sample[p.mum.id].data["GQ"] if p.mum and p.mum.id in rec.call_for_sample else None)
                        gq_d = getMax(rec.call_for_sample[p.dad.id].data["GQ"] if p.dad and p.dad.id in rec.call_for_sample else None)
                        if gq_m is not None and gq_m > min_gq_dnm and gq_d is not None and gq_d > min_gq_dnm and gq > min_gq_dnm:
                            if gt_m == 0 and gt_d == 0:
                                # get additional metadata for DNM filtering
                                dp = call.data["DP"] if "DP" in call.data else None
                                ad = call.data["AD"] if "AD" in call.data else None
                                ab = getMax(call.data["AB"]) if "AB" in call.data else None
                                aq = getMax(rec.INFO["AQ"]) if "AQ" in rec.INFO else None
                                filtered=False
                                # check for minimum AD of alt-allele if AD field set
                                if ad is not None and len(ad)>0 and ad[1]<4:
                                    filtered=True
                                # check for allelic balance
                                if ab is not None and ( ab < 0.2 or ab > 0.8):
                                    filtered=True
                                # check for allele quality
                                if aq is not None and aq < 30:
                                    filtered=True
                                if not filtered:
                                    sup_dnm += 1 # supports DNM
                            else:
                                sup_dnm -= 1 # inherited
                                
                            if ( gt_m > 0 and p.mum not in self.affected_ids() ) or ( gt_d > 0 and p.dad not in self.affected_ids() ):
                                sup_dom -= 1
                            elif ( gt_m > 0 and p.mum in self.affected_ids() ) or ( gt_d > 0 and p.dad in self.affected_ids() ):
                                sup_dom +=1
                elif gt==2: # hom
                    sup_dom -= 1
                    if is_inherited:
                        sup_dnm -= 1 
                    if ( gt_m is None or gt_m ==1 ) and ( gt_d is None or gt_d ==1 ): # hom inherited from mum&dad
                        sup_rec += 1
                    else:
                        sup_rec -= 1
            else:
                if gt==0: # homref
                    pass
                elif gt==1: # het
                    sup_dnm -= 1
                    sup_dom -= 1
                elif gt==2: # hom
                    sup_dnm -= 1
                    sup_rec -= 1
                    sup_dom -= 1
        gt_str = "\t".join(genotypes)
        gq_str = "\t".join(genotype_quals)
        return [sup_rec, sup_dom, sup_dnm, num_hq_calls,hom_aff,hom_unaff,het_aff,het_unaff, gt_str, gq_str, sup_comhet_mum, sup_comhet_dad]


    def __repr__(self):
        return self.__str__()  

    def __str__(self):
        return ('ped %s [%s]' % (self.name, ", ".join(str(x) for x in self.families.values())))

    
def add_link(ped2a2b2ids, ped_name, a, b, id):
    if ped_name not in ped2a2b2ids:
        ped2a2b2ids[ped_name]={}
    if a not in ped2a2b2ids[ped_name]:
        ped2a2b2ids[ped_name][a]={}
    if b not in ped2a2b2ids[ped_name][a]:
        ped2a2b2ids[ped_name][a][b]=[]
    ped2a2b2ids[ped_name][a][b].append(id)
    return ped2a2b2ids

def add_link2(ped2gid2a2b2ids, ped_name, gid, a, b, id):
    if ped_name not in ped2gid2a2b2ids:
        ped2gid2a2b2ids[ped_name]={}
    if gid not in ped2gid2a2b2ids[ped_name]:
        ped2gid2a2b2ids[ped_name][gid]={}
    if a not in ped2gid2a2b2ids[ped_name][gid]:
        ped2gid2a2b2ids[ped_name][gid][a]={}
    if b not in ped2gid2a2b2ids[ped_name][gid][a]:
        ped2gid2a2b2ids[ped_name][gid][a][b]=[]
    ped2gid2a2b2ids[ped_name][gid][a][b].append(id)
    return ped2gid2a2b2ids

# calculates what fraction of a is covered by b. Works only for single pyranges intervals!
def calc_region_coverage(a, b):
    #print(a)
    #print(b)
    if len(a.index)!=1 or len(b.index)!=1:
        raise Exception("multiple pr objects passed")
    if a.Chromosome[0] != b.Chromosome[0]:
        return 0.0
    len_a=a.End-a.Start+1
    d=a.subtract(b)
    if len(d)==0:
        return 1.0
    len_d=d.End-d.Start+1
    return (1-len_d/len_a)

def get_output_rec(inrec, include_info_fields, filter, out_id, info_keep=["SVTYPE", "SVLEN", "POS", "END"], info_add=OrderedDict()):
    rec=copy.deepcopy(inrec)
    if filter:
        rec.add_filter(filter)
    if not include_info_fields:
        for f in info_keep:
            if f in rec.INFO:
                info_add[f]=rec.INFO[f]
        rec.INFO.clear()
    if out_id:
        rec.INFO["var2reg_id"]=out_id
    if info_add:
        rec.INFO.update(info_add)
    return rec

#============================================================================
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-conf","--config", type=existing_file, required=True, dest="confF", help="JSON config file")
parser.add_argument("-t","--threads", type=int, required=False, default=4, dest="threads",  help="Number of threads (default: 4)")
parser.add_argument("-o","--out", type=str, required=False, dest="outdir", metavar="outdir", help="Output folder")
args = parser.parse_args() 


print(logo)
startTime = time.time()

# ignore warnings (mainly from vcfpy)
if not sys.warnoptions:
    warnings.filterwarnings(action="ignore", module="vcfpy")
#============================================================================

config_file = args.confF
outdir = args.outdir if args.outdir else os.getcwd()
threads = args.threads

# load + check config
config = json.load(open(config_file), object_pairs_hook=OrderedDict)
check_config(config)

# get chunk size
chunksize = config["input_data"]["chunksize"] if "chunksize" in config["input_data"] else 1000000
print("Setting chunksize to %i" % (chunksize))

# ensure dirs
if not outdir.endswith("/"):
    outdir += "/"
outdir += config["dataset_name"] + "/"
if not os.path.exists(outdir):
    print("Creating dir " + outdir)
    os.makedirs(outdir)
log = os.path.join(outdir, pipename.replace(" ", "-") + ".log")

# debug mode?
debug = False
if "debug" in config and config["debug"]:
    debug = True

# add commandline params to config and write to output dir
config["CMD"]={}
for arg in vars(args):
    config["CMD"][arg] = getattr(args, arg) 
# write effective config
outF=outdir + config["dataset_name"] + ".v2r.effective_conf.json"
with open(outF, 'w') as outs:
    print(json.dumps(config, indent=4, sort_keys=False), file=outs)

chroms = CHR_hg38 = ['chr'+str(x) for x in range(1,23)]+["chrX", "chrY"]
if "chrom" in config:
    chroms = config["filters"]["chrom"]
print("Handling chromosomes %s" % (", ".join(chroms)))

# read gene name alias table
alias2gid={}
if "alias_table" in config["input_data"]:
    d = pd.read_csv(config["input_data"]["alias_table"],delimiter='\t',encoding='utf-8')
    alias2gid=pd.Series(d.symbol.values,index=d.alias).to_dict()
print("Loaded %i gene symbol aliases" % (len(alias2gid)))

# check model output
max_var_per_cat = config["input_data"]["max_var_per_cat"] if "max_var_per_cat" in config["input_data"] else None
include_na_models = config["input_data"]["include_na_models"] if "include_na_models" in config["input_data"] else True
if not include_na_models:
    print("===========================================================================")
    print("NOTE that NA model entries will be filtered with the current configuration!")
    print("===========================================================================")

# read SO terms and weights
print("Reading SO terms + weights")
terms={}
for t in config["so_term"]:
    d_score=config["so_term"][t]["d_score"]
    a_score=config["so_term"][t]["a_score"]
    a_type=config["so_term"][t]["a_type"]
    so_id=config["so_term"][t]["so_id"]
    descr=config["so_term"][t]["descr"]
    terms[t] = SoTerm(t,d_score,a_score,a_type,so_id,descr)
print("\t%s" % (terms))

# read gene annotations
gene_anno={}
gene_anno_columns=None
if "gene_anno_table" in config["input_data"]:
    nline = sum(1 for i in gzip.open(config["input_data"]["gene_anno_table"], 'rb'))
    print("Accessing file with %i additional gene annotations" % (nline))
    with tqdm(total=nline, unit=' annotations', desc='Loading gene annotations') as bar:  
        for d in pd.read_csv(config["input_data"]["gene_anno_table"],
                                 delimiter='\t',
                                 encoding='utf-8',
                                 comment="#",
                                 float_precision='high',
                                 chunksize=chunksize,
                                 quoting=csv.QUOTE_NONE,
                                 index_col=False,
                                 error_bad_lines=False): 
            gene_anno_columns=d.columns.tolist()
            del gene_anno_columns[0]
            for index, row in d.iterrows():
                l = row.tolist()
                del l[0]
                l = ["NA" if math.isnan(x) else x for x in l]
                gene_anno[row[0]]=l
                bar.update(len(gene_anno))
        bar.update(nline)
    print("\tread %i gene annotations" % (len(gene_anno)))        
else:
    print("No gene annotations configured!")


# =================================================================================
# read all regulatory regions
# =================================================================================
nline = sum(1 for i in gzip.open(config["input_data"]["reg_db"], 'rb'))
print("Accessing file with %i regulatory regions" % (nline))
with tqdm(total=nline, unit=' reg_regions', desc='Loading regulatory regions annotations') as bar:  
    names = ['Chromosome', 'Start', 'End', 'id', 'std_type', 'db_source', 'closestGene_symbol', 'closestGene_dist', 'affected_genes', 'N_methods'] 
    dtype = {'Chromosome': str, 'Start': int, 'End': int, 'id': str, 'std_type': str, 'db_source': str, 'closestGene_symbol': str, 'closestGene_dist': int, 'affected_genes': str, 'N_methods': str}
    reg_regions=[]
    id2reg_regions={}
    for d in pd.read_csv(config["input_data"]["reg_db"],
                             delimiter='\t',
                             encoding='utf-8',
                             header=None,
                             names=names,
                             dtype=dtype,
                             comment="#",
                             float_precision='high',
                             chunksize=chunksize,
                             quoting=csv.QUOTE_NONE,
                             index_col=False,
                             low_memory=False,
                             error_bad_lines=False,
                             keep_default_na=False): # this is required to interpret empty strings as NA instead of NaN 
        if chroms:
            d = d[d['Chromosome'].isin(chroms)]
        d['Start']+=1 # 0-based coords!
        reg_regions.append(pr.PyRanges(d))
        bar.update(len(d.index))
    bar.update(nline)
print("\tbuilding index") 
reg_regions=pr.concat(reg_regions)
for index, region in reg_regions.df.iterrows():
    id2reg_regions[region['id']]=region 
print("\tread %i regulatory regions" % (len(reg_regions)))        
with pd.option_context('display.max_rows',5):
    print(reg_regions.df)
    
# =================================================================================
# read gene annotations from GFF
# =================================================================================
print("Reading gene annotations")
gff = pr.read_gff3(config["input_data"]["gene_gff"])
exon_int = gff[gff.Feature == 'exon']
utr3_int = gff[gff.Feature == 'three_prime_UTR']
utr5_int = gff[gff.Feature == 'five_prime_UTR']
utr_int = pr.concat([utr3_int, utr5_int])
gene_int = gff[gff.Feature == 'gene']
print("Loaded %i exon, %i utr and %i gene annotations from %s" % (len(exon_int), len(utr_int), len(gene_int), config["input_data"]["gene_gff"]))

# read samples we have data for from VCF header
reader = vcfpy.Reader.from_path(config["input_data"]["snv_vcf"])
reader2 = vcfpy.Reader.from_path(config["input_data"]["cnv_vcf"]) if config["input_data"]["cnv_vcf"] else None


# =================================================================================
# read input pedigrees
# =================================================================================
print("Reading pedigrees and gene rank data")
peds=[]
gados=[]
exos=[]
for id in config["input_data"]["pedigrees"]:
    pedf = config["input_data"]["ped_pattern"].format(id=id)
    gadof = config["input_data"]["gado_pattern"].format(id=id)
    anyexo = False
    for im in EXOMISER_IM:
        exof = config["input_data"]["exomiser_pattern"].format(IM=im, id=id)
        anyexo = anyexo or files_exist(exof)
    if files_exist([pedf, gadof]) and anyexo:
        ped = Pedigree(pedf)
        # do we have data for at least one sample?
        if len(list(set(ped.all_ids()) & set(reader.header.samples.names))) > 0:
            peds += [ped]
            gid2rnk={}
            d = pd.read_csv(gadof, delimiter='\t', encoding='utf-8')
            for index, row in d.iterrows():
                gid = row['Hgnc']
                gid = alias2gid[gid] if gid in alias2gid else gid
                gid2rnk[gid]=row['Zscore']
            gados+=[gid2rnk]
            exo2rnk={}
            # read exomiser tables
            for im in EXOMISER_IM:
                exof = config["input_data"]["exomiser_pattern"].format(IM=im, id=id)
                if files_exist([exof]):
                    d = pd.read_csv(exof, delimiter='\t', encoding='utf-8')
                    for index, row in d.iterrows():
                        exo2rnk[row['#GENE_SYMBOL']]=row['EXOMISER_GENE_PHENO_SCORE']
            exos+=[exo2rnk]       
            print("\t%s: %i gado_rnk, %i exo_ps" % (ped.name, len(gid2rnk), len(exo2rnk)))
        else:
            print("\tDid not find any data for the following samples in the snv VCF: %s. Ignoring ped %s." % (",".join(ped.all_ids()), id))
    else:
        print("\tCannot load all required data for id %s. Some files (ped, gado or exomiser results) do not exist. Ignoring (%s, %s, %s)" % (id, pedf, gadof, config["input_data"]["exomiser_pattern"]) )
    #print("\tsamples: %s" % (ped.all_ids()))
    #print("\taffected_samples: %s" % (ped.affected_ids()))
print("\tDone. Read %i pedigrees and %i gado rank lists" % (len(peds), len(gados)))
if len(peds) == 0:
    sys.exit("Check input configuration / data...")


# =================================================================================
# split VCFs
# =================================================================================
print("Splitting input VCFs per chrom")

def split_vcf_by_chrom(vcf_file, chroms, out_prefix, chunksize=10000, threads=1):
    # read VCF header via vcfpy (# full file: 13Mio SNVs)
    reader = vcfpy.Reader.from_path(vcf_file)
    fout={}
    total_var_count=0 # for counting # of vars
    for c in chroms:
        fout[c] = out_prefix+"."+str(c)+".vcf"
        writer = vcfpy.Writer.from_path(fout[c], reader.header)
        writer.close()
    try:
        for d in pd.read_csv(vcf_file,
                             delimiter='\t',
                             encoding='utf-8',
                             header=None,
                             dtype=str,
                             comment="#",
                             float_precision='high',
                             chunksize=chunksize,
                             quoting=csv.QUOTE_NONE,
                             error_bad_lines=False): 
            found_chr = set(d[0])
            for c in found_chr:
                if c in fout:
                    total_var_count+=len(d.index)
                    fil = d[d.iloc[:,0]==c]
                    fil.to_csv(fout[c], sep='\t', mode='a', header=None, index=False)
        for c in chroms:
            bgzip(fout[c], index=True, delinFile=True, maxthreads=threads)
            fout[c] = fout[c]+".gz"
    except pd.io.common.EmptyDataError:
        print("No data in VCF file %s" % vcf_file)
    return(fout, total_var_count)

snvF_chr, total_var_count_snv = split_vcf_by_chrom(config["input_data"]["snv_vcf"], chroms, outdir+"snv", chunksize, threads)
cnvF_chr, total_var_count_cnv = split_vcf_by_chrom(config["input_data"]["cnv_vcf"], chroms, outdir+"cnv", chunksize, threads) if config["input_data"]["cnv_vcf"] else None

# =================================================================================
# parallel iterate variants per chrom and write var2reg tables
# =================================================================================
def process_chrom(args):
    
    chrom,  snvF, cnvF, peds, gados, exos, gene_anno, gene_anno_columns, max_var_per_cat, include_na_models, chunksize, debug, outdir, queue  = args 
    tqdm.write("Parsing %s SNVs from %s and CNVs from %s" % (chrom, snvF, cnvF))
     
    # for storing vars per inheritance model
    ped2gid2mod2vid={}
    # for storing comphet candidates
    ped2gid2vid2mumordad2pid={}
    # SO terms that were ignored
    ignored_SO_terms = set()
    # for keeping track of current variant
    wrote2results={} 
    # variant and record ids
    var_id_counter=0
    rec_id_counter=0

    # output records
    records={}
    for ped in peds:
        records[ped.name] = list()

    # output files
    recF=OrderedDict()
    recOut={}
    for ped in peds:
        recF[ped.name] = outdir+"."+chrom+"." + ped.name + ".rec_tmp"
        recOut[ped.name] = open(recF[ped.name], 'w')
    logF=outdir+"."+chrom+".log_tmp"
    
    # filter counts
    filter_counts={}
    total_var_count={}
    for vcf_type in ["SNV", "CNV"]:
        filter_counts[vcf_type]={}
        total_var_count[vcf_type]=0
        for tag in filter_tags:
            filter_counts[vcf_type][tag]=0
    
    # iterate SNVs and CNVs
    with open(logF, 'w') as log:
        for vcf_type in ["SNV", "CNV"]:
            #print("Analyzing %s records" % (vcf_type))
            # for warning about slow pyranges queries
            did_warn_pyranges=False
        
            # load data
            infile = snvF if vcf_type == "SNV" else cnvF
            if not infile:
                tqdm.write("No VCF file configured for type %s." % infile )
                continue # no CNV file configured
            
            # read VCF header via vcfpy (# full file: 13Mio SNVs)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore") # supress warnings from incomplete headers
                header = reader.header if vcf_type=="SNV" else reader2.header
                parser = vcfpy.parser.RecordParser(header, header.samples)
                if vcf_type=="CNV":
                        # fix issue with missing CN header in CNV vcf
                        header.add_format_line(vcfpy.OrderedDict([('ID', 'CN'), ('Description', 'Copynumber'), ('Number', '1'), ('Type', 'Float')]))
    
            # iterate variants
            names = ['CHROM', 'POS', 'ID', 'REF', 'ALT',
                       'QUAL', 'FILTER', 'INFO', 'FORMAT'] 
            dtype = {'CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT':str}
            for s in header.samples.names:
                names = names + [s]
                dtype[s] = str
        
            var_calls = 0
            for d in pd.read_csv(infile,
                                 delimiter='\t',
                                 encoding='utf-8',
                                 header=None,
                                 names=names,
                                 dtype=dtype,
                                 comment="#",
                                 float_precision='high',
                                 chunksize=chunksize,
                                 quoting=csv.QUOTE_NONE,
                                 error_bad_lines=False): 
                # max var per cat
                if max_var_per_cat and var_calls >= max_var_per_cat:
                    break
            
                # iterate only pass variants on given chrom
                all = len(d.index)
                total_var_count[vcf_type]+=all
                d = d.query("(FILTER=='PASS' | FILTER=='.') & CHROM=='"+chrom+"'")
                filter_counts[vcf_type]["NP"]+=all-len(d.index)
                 
                for index, row in d.iterrows():
                    var_id_counter += 1
                    var_calls += 1
    
                    # max var per cat
                    if max_var_per_cat and var_calls >= max_var_per_cat:
                        break
            
                    # parse VCF record
                    rec = parser.parse_line("\t".join(map(str, row)))
            
                    # get pos/len
                    altstr=','.join(map(serialize, rec.ALT))
                    alt1=altstr.split(',')[0]
                    pos, end=None, None
                    if vcf_type == "SNV":
                        pos = rec.POS
                        end = pos + len(alt1) - 1 
                    elif vcf_type=="CNV":
                        if ("POS" in rec.INFO) and ("END" in rec.INFO):
                            pos = int(rec.INFO.get("POS"))
                            end = int(rec.INFO.get("END"))
                    if not pos or not end:
                        # INCOMPLETE_DATA. ignore.
                        if debug:
                            print("ID\t%s:%s" % (rec.CHROM, rec.POS), file=log)
                        filter_counts[vcf_type]["ID"]+=1
                        continue
                        
                    # get vartype
                    vartype = None
                    if vcf_type == "SNV":
                        if len(rec.REF)==len(alt1):
                            vartype = "SNV" if len(rec.REF)==1 else "MNV"
                        else:
                            vartype="INDEL"
                    elif vcf_type=="CNV":
                        # calculate variant type
                        is_unspecific=False
                        for a in rec.ALT:
                            if a.type == "BND":  # skip general breakpoint annotations as too unspecific 
                                is_unspecific = True
                            vartype = a.value if a.value != "None" else a.type
                        if vartype is None or is_unspecific:
                            # UNSPECIFIC_TYPE. ignore
                            if debug:
                                print("UT\t%s:%s" % (rec.CHROM, rec.POS), file=log)
                            filter_counts[vcf_type]["UT"]+=1
                            continue
                        # this will map, e.g., DEL=>DEL, DEL:ME=>DEL:ME, DEL:ME:something=>DEL:ME
                        vartype = ":".join(vartype.split(":")[0:2])
                        
                    # filter by max pop AF
                    #maxAF = calc_max_AF(rec, ["global_A1000G_AF", "global_ExAC_AF", "global_UK10K_AF", "global_gnomAD_AF"] if vcf_type =="SNV" else ["GD_POPMAX_AF", "1000g_max_AF"])
                    maxAF = calc_score(rec, config, "max_pop_af_snv", 0) if vcf_type =="SNV" else calc_score(rec, config,"max_pop_af_cnv",0)
                    if maxAF is not None and maxAF > config["filters"]["max_pop_af"]:
                        # TOO_COMMON. ignore
                        if debug:
                            print("AF\t%s:%s" % (rec.CHROM, rec.POS), file=log)
                        filter_counts[vcf_type]["AF"]+=1
                        continue
        
                    # cohort freq: FIXME (better impl for CNVS based on partial overlap?)
                    cohort_calls = 0
                    all_calls = 0
                    for call in rec.calls:
                        gt = call.gt_type # hom_ref = 0, het = 1, hom_alt = 2 (which alt is untracked), uncalled = None
                        gq = call.data["GQ"]
                        maxgq = getMax(gq)
                        if maxgq is None:
                            maxgq = 0    
                        if gt is not None and maxgq > config["filters"]["min_gq"]:
                            cohort_calls += gt
                            all_calls += 2
                    cohortAF = 0 if all_calls == 0 else float(cohort_calls) / float(all_calls) # allele frequency in cohort
                    
                    # iterate over all pedigrees
                    is_nocall = 0
                    is_ni = 0
                    is_low_score = 0
                    for idx, ped in enumerate(peds):
                        wrote2results[ped.name]=False
                    
                        # calc support for inh models
                        min_gq = config["filters"]["min_gq"] if vcf_type == "SNV" else -1 # GQ is not a reliable measure for CNVs
                        min_gq_dnm = (config["filters"]["min_gq_dnm"] if "min_gq_dnm" in config["filters"] else 30) if vcf_type == "SNV" else -1 # GQ is not a reliable measure for CNVs
                        sup_rec, sup_dom, sup_dnm, num_hq_calls, hom_aff, hom_unaff, het_aff, het_unaff, gt_str, gq_str, sup_comhet_mum, sup_comhet_dad = ped.calc_inh_support(rec, min_gq, min_gq_dnm )
                        max_sup = getMax([sup_rec, sup_dom, sup_dnm])
                        sup_comhet = len(sup_comhet_mum) + len(sup_comhet_dad)
                        num_affected_calls = hom_aff + het_aff
                        
                        # no call with GQ>thresh in affected samples. ignore.
                        if num_affected_calls == 0: 
                            is_nocall+=1
                            continue
                        # no support for any inh model or comphet. ignore.
                        if not include_na_models and sup_comhet == 0 and ( max_sup is None or max_sup <= 0):
                            is_ni+=1
                            continue
                        # collect INFO values   
                        infos=[]
                        for f in config["output_fields"]["included_info_fields"]:
                            infos.append(fmt(getMax(str(rec.INFO.get(f, "NA")).split(","))))
                        for f in config["output_fields"]["included_info_factors"]:
                            infos.append(str(rec.INFO.get(f, "NA")))
                        for f in config["output_fields"]["included_info_flags"]:
                            infos.append("1" if f in rec.INFO else "0")
                            
                        #--------------------------------
                        # Find overlap with genic regions
                        #--------------------------------
                        # get term with max weight per gene
                        gene2term={}
                        gene2HGVSp={} 
                        if vcf_type == "SNV":
                            if "ANN" in rec.INFO:
                                for a in rec.INFO["ANN"]:
                                    # 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO
                                    # ex: ANN=T|splice_donor_variant&intron_variant|HIGH|HDLBP|ENSG00000115677|transcript|ENST00000310931.8|protein_coding|14/27|c.1731+1G>A||||||,T|....
                                    for aa in rec.ALT:
                                        dat = a.split("|")
                                        gid=dat[3] # gene name
                                        gid = alias2gid[gid] if gid in alias2gid else gid
                                        #print(dat[0]+"/"+aa.serialize() + "-" + dat[1])
                                        if dat[0] == aa.serialize(): # check right alt-allele
                                            for imp in dat[1].split("&"): # iterate annotations
                                                if imp in terms:
                                                    if gid in gene2term:
                                                        # replace if more deleterious only (i.e., respect input order)
                                                        if terms[gene2term[gid]].d_score<terms[imp].d_score:
                                                            #print("%s < %s" %(gene2term[gid], imp))
                                                            gene2term[gid] = imp
                                                        #else:
                                                        #    print("%s >= %s" %(gene2term[gid], imp))
                                                    else:
                                                        gene2term[gid] = imp
                                                else:
                                                    ignored_SO_terms.add(imp)
                                            if dat[10]: # iterate HGVS.p
                                                hgvsp = gene2HGVSp[gid] if gid in gene2HGVSp else []
                                                hgvsp += [ (dat[6]+":" + dat[10]).replace(',', '_') ]
                                                gene2HGVSp[gid] = hgvsp
                            else:
                                warnings.warn("No ANN field found in SNV entry at %s:%i" % (rec.CHROM, pos) )
                        elif vcf_type == "CNV":
                                query = pr.PyRanges(chromosomes=[rec.CHROM], starts=[min(pos, end)], ends=[max(pos, end)])
                                join = query.join(exon_int)
                                # FIXME: select most severe per gene!
                                if 'gene_name' in join.columns:
                                    overlapping_exons = set(join.gene_name.values)   
                                    for gid in overlapping_exons:
                                        gene2term[gid]="exonic_sv" # must be in config file!
                                join = query.join(utr_int)
                                if 'gene_name' in join.columns:
                                    overlapping_utrs = set(join.gene_name.values)   
                                    for gid in overlapping_utrs:
                                        if gid not in gene2term:
                                            gene2term[gid]="utr_sv"
                                join = query.join(gene_int)
                                if 'gene_name' in join.columns:
                                    overlapping_introns = set(join.gene_name.values)   
                                    for gid in overlapping_introns:
                                        if gid not in gene2term:
                                            gene2term[gid]="intronic_sv" # we assume 'intron' here if not overlapping with exon/utr annotations of this gene.                     
        
                        handled_gids=[]
                        for gid in gene2term:
                            imp = gene2term[gid] # impact
                            term = terms[gene2term[gid]] # SO:term
                            aa_change = gene2HGVSp[gid] if gid in gene2HGVSp else ["NA"]
                            d_score = term.d_score # default
                            a_score = term.a_score
                            # FIXME: how to calc these scores for SVs?
                            if vcf_type == "SNV":
                                if term.type == "intron":
                                    # intronic variant, calculate d_score from NC scores or spliceAI et al.
                                    d_score=calc_score(rec, config, "intron", d_score)
                                elif term.type == "splicing":
                                    # splicing region variant according to VEP (3bp in exon, 8bp in intron). calc weight from spliceAI et al.
                                    d_score=calc_score(rec, config, "splicing", d_score)
                                elif term.type == "utr":
                                    # UTR variant. Calculate from NC scores
                                    d_score=calc_score(rec, config, "utr", d_score)
                            
                            if vcf_type == "CNV" or d_score >= config["filters"]["min_dscore"]:
                                handled_gids += gid
                                
                                rec_id = chrom + ":" + str(rec_id_counter)
                                var_id = chrom + ":" + str(var_id_counter)
                                records[ped.name].append(  [rec_id, var_id, gid, rec.CHROM, pos, end, rec.REF, 
                                    altstr, vartype, ",".join(aa_change), ",".join(rec.ID) if len(rec.ID)>0 else "NA", fmt(maxAF), fmt(cohortAF), 
                                    "NA", "anno", term.type, "NA", term.term, fmt(d_score), fmt(a_score), 
                                    sup_rec, sup_dom, sup_dnm, 
                                    ",".join(sup_comhet_mum) if sup_comhet_mum else "NA", ",".join(sup_comhet_dad) if sup_comhet_dad else "NA",  
                                    num_hq_calls, hom_aff, hom_unaff, het_aff, het_unaff,
                                    gt_str, gq_str, "\t".join(infos) ] )
                                rec_id_counter += 1
                                
                                # add gene<>inh_model<>var links
                                if sup_rec > 0:
                                    ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "recessive", rec_id)
                                if sup_dom > 0:
                                    ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "dominant", rec_id)
                                if sup_dnm > 0:
                                    ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "denovo", rec_id)
                                if include_na_models:
                                    ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "NA", rec_id)
            
                                # comphet. sup_comhet_mum, sup_comhet_dad contain the pids of affected samples that inherit a het from mum/dad
                                for pid in sup_comhet_mum:
                                    ped2gid2vid2mumordad2pid = add_link2(ped2gid2vid2mumordad2pid, ped.name, gid, rec_id, "mum", pid)
                                for pid in sup_comhet_dad:
                                    ped2gid2vid2mumordad2pid = add_link2(ped2gid2vid2mumordad2pid, ped.name, gid, rec_id, "dad", pid)
            
                                wrote2results[ped.name]=True
                            
                        #--------------------------------
                        # Find overlap with regulatory regions. Filter for dscore only for SNVs
                        #--------------------------------    
                        if vcf_type == "SNV":
                            d_score=calc_score(rec, config, "noncoding", 0)
                        else:
                            d_score=None
                        if vcf_type == "CNV" or d_score >= config["filters"]["min_dscore"]:
                            gid2reg={}
                            gid2assoctype={}
                            if vcf_type == "SNV" and "snv_reg_id_info_field" in config["filters"]:
                                # try to get reg_id values from INFO fields
                                if rec.INFO.get(config["filters"]["snv_reg_id_info_field"], None):
                                    for rid in rec.INFO.get(config["filters"]["snv_reg_id_info_field"]).split(","):
                                        if rid in id2reg_regions:
                                            region=id2reg_regions[rid]
                                            has_assoc = region['affected_genes'] and len(str(region['affected_genes']))>0
                                            gids = str(region['affected_genes']).split(",") if has_assoc else str(region['closestGene_symbol']).split(",")
                                            assoc_type = "reg_db" if has_assoc else "closest_gene"
                                            for gid in gids:
                                                gid = alias2gid[gid] if gid in alias2gid else gid
                                                if gid not in gid2reg:
                                                    gid2reg[gid]=[]
                                                gid2reg[gid].append(region)
                                                gid2assoctype[gid] = assoc_type
                            elif vcf_type == "CNV" and "cnv_reg_id_info_field" in config["filters"]:
                                # try to get reg_id values from INFO fields
                                if rec.INFO.get(config["filters"]["cnv_reg_id_info_field"], None):
                                    for rid in rec.INFO.get(config["filters"]["cnv_reg_id_info_field"]).split(","):
                                        if rid in id2reg_regions:
                                            region=id2reg_regions[rid]
                                            has_assoc = region['affected_genes'] and len(str(region['affected_genes']))>0
                                            assoc_type = "reg_db" if has_assoc else "closest_gene"
                                            gids = str(region['affected_genes']).split(",") if has_assoc else str(region['closestGene_symbol']).split(",")
                                            for gid in gids:
                                                gid = alias2gid[gid] if gid in alias2gid else gid
                                                if gid not in gid2reg:
                                                    gid2reg[gid]=[]
                                                gid2reg[gid].append(region)
                                                gid2assoctype[gid] = assoc_type
                            else:
                                # pyranges query
                                if not did_warn_pyranges:
                                     warnings.warn("Could not extract pre-annotated reg ids from VCF. Searching for regulatory regions via pyranges which is much slower!")
                                     did_warn_pyranges=True
                                query=pr.PyRanges(chromosomes=rec.CHROM, starts=[pos], ends=[end])
                                overlapping_reg_regions = reg_regions.join(query)
                                df=overlapping_reg_regions.df
                                if df is not None:
                                    for index, region in df.iterrows():
                                        has_assoc = region['affected_genes'] and len(str(region['affected_genes']))>0
                                        assoc_type = "reg_db" if has_assoc else "closest_gene"
                                        gids = str(region['affected_genes']).split(",") if has_assoc else str(region['closestGene_symbol']).split(",")
                                        for gid in gids:
                                            gid = alias2gid[gid] if gid in alias2gid else gid
                                            if gid not in gid2reg:
                                                gid2reg[gid]=[]
                                            gid2reg[gid].append(region)
                                            gid2assoctype[gid] = assoc_type
                            #--------------------------------
                            # write one var per gene/region
                            #--------------------------------    
                            for gid in gid2reg.keys():
                                for region in gid2reg[gid]: 
                                    a_score = 0
                                    if region['db_source'] in config["def_a_score"]:
                                        # read default a_score from config
                                        a_score = config["def_a_score"][region['db_source']]
                                    rec_id = chrom + ":" + str(rec_id_counter)
                                    var_id = chrom + ":" + str(var_id_counter)
                                    records[ped.name].append(  
                                        [rec_id, var_id, gid, rec.CHROM, pos, end, rec.REF, 
                                         altstr, vartype, "NA",",".join(rec.ID) if len(rec.ID)>0 else "NA", fmt(maxAF), fmt(cohortAF), 
                                         region['id'], region['db_source'], gid2assoctype[gid], region['closestGene_dist'], region['std_type']+"_variant", 
                                         fmt(d_score), fmt(a_score), 
                                         sup_rec, sup_dom, sup_dnm,
                                         ",".join(sup_comhet_mum) if sup_comhet_mum else "NA", ",".join(sup_comhet_dad) if sup_comhet_dad else "NA",
                                         num_hq_calls, hom_aff, hom_unaff, het_aff, het_unaff,
                                         gt_str, gq_str, "\t".join(infos)] )
                                    rec_id_counter += 1
                                        
                                    # add gene<>inh_model<>var links
                                    if sup_rec > 0:
                                        ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "recessive", rec_id)
                                    if sup_dom > 0:
                                        ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "dominant", rec_id)
                                    if sup_dnm > 0:
                                        ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "denovo", rec_id)
                                    if include_na_models:
                                        ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "NA", rec_id)
            
                                    # comphet. sup_comhet_mum, sup_comhet_dad contain the pids of affected samples that inherit a het from mum/dad
                                    for pid in sup_comhet_mum:
                                        ped2gid2vid2mumordad2pid = add_link2(ped2gid2vid2mumordad2pid, ped.name, gid, rec_id, "mum", pid)
                                    for pid in sup_comhet_dad:
                                        ped2gid2vid2mumordad2pid = add_link2(ped2gid2vid2mumordad2pid, ped.name, gid, rec_id, "dad", pid)
            
                                    wrote2results[ped.name]=True
                        # LOW_SCORE variants. ignore
                        if not wrote2results[ped.name]:
                            is_low_score+=1
                            
                    if is_nocall == len(peds):
                        # NO_CALL. ignore
                        if debug:
                            print("NC\t%s:%s" % (rec.CHROM, rec.POS), file=log)
                        filter_counts[vcf_type]["NC"]+=1
                    if is_ni == len(peds):
                        # NO_INH_MODEL. ignore
                        if debug:
                            print("NI\t%s:%s" % (rec.CHROM, rec.POS), file=log)
                        filter_counts[vcf_type]["NI"]+=1
                    if is_low_score == len(peds):
                        if debug:
                            print("LS\t%s:%s" % (rec.CHROM, rec.POS), file=log)
                        filter_counts[vcf_type]["LS"]+=1
                    
                    # progress bar
                    queue.put(1)
                           
                # write records to disc
                if len(records) > 100:
                    for ped_id in records.keys():
                        for rdata in records[ped_id]:
                            print( "\t".join( str(x) for x in rdata), file=recOut[ped_id])
                        records[ped_id]=list()
                          
        # write remaining records       
        for ped_id in records.keys():
            for rdata in records[ped_id]:
                print( "\t".join( str(x) for x in rdata), file=recOut[ped_id])
            records[ped_id]=list()    
                
    # close streams
    for ped in peds:
        recOut[ped.name].close()
    #print("Done iterating %s" % chrom) 
            
    # write gene tables
    gene_ids = {} 
    records_genes = {}
    gene_count=0
    for idx, ped in enumerate(peds):
        records_genes[ped.name] = list()
        gene_ids[ped.name]=0
        if ped.name in ped2gid2mod2vid:
            for gid in ped2gid2mod2vid[ped.name]:
                for ms in ped2gid2mod2vid[ped.name][gid]:
                    gado_z = str(gados[idx][gid]) if gid in gados[idx] else "NA"
                    exo_ps = str(exos[idx][gid]) if gid in exos[idx] else "NA"
                    annostr=""
                    if gene_anno_columns:
                        anno=gene_anno[gid] if gid in gene_anno else ["NA"] * len(gene_anno_columns)
                        annostr="\t".join(str(x) for x in anno)
                    records_genes[ped.name].append([gene_ids[ped.name], gid, gado_z, exo_ps, ms, len(ped2gid2mod2vid[ped.name][gid][ms]), ped2gid2mod2vid[ped.name][gid][ms], annostr])
                    gene_count+=1
                gene_ids[ped.name] += 1

    # write comphet table
    comphet_ids = {} 
    records_comphet = {}
    comphet_count=0
    for idx, ped in enumerate(peds):
        records_comphet[ped.name] = list()
        comphet_ids[ped.name]=0
        if ped.name in ped2gid2vid2mumordad2pid:
            for gid in ped2gid2vid2mumordad2pid[ped.name]:
                vids = ped2gid2vid2mumordad2pid[ped.name][gid]
                for i1, v1 in enumerate(vids):
                    for i2, v2 in enumerate(vids):
                        if v1 < v2:
                            # samples that inherit v1 from mum 
                            mum1 = ped2gid2vid2mumordad2pid[ped.name][gid][v1]["mum"] if "mum" in ped2gid2vid2mumordad2pid[ped.name][gid][v1] else []
                            dad1 = ped2gid2vid2mumordad2pid[ped.name][gid][v1]["dad"] if "dad" in ped2gid2vid2mumordad2pid[ped.name][gid][v1] else []  
                            mum2 = ped2gid2vid2mumordad2pid[ped.name][gid][v2]["mum"] if "mum" in ped2gid2vid2mumordad2pid[ped.name][gid][v2] else []
                            dad2 = ped2gid2vid2mumordad2pid[ped.name][gid][v2]["dad"] if "dad" in ped2gid2vid2mumordad2pid[ped.name][gid][v2] else []
                            m1d2 = list(set(mum1) & set(dad2)) # samples for which there is a hit inherited from mum and dad
                            m2d1 = list(set(mum2) & set(dad1))  
                            candidate_pid = set(m1d2 + m2d1)
                            if candidate_pid:
                                records_comphet[ped.name].append([comphet_ids[ped.name], 
                                                          gid, 
                                                          v1,
                                                          v2,
                                                          len(candidate_pid), 
                                                          ",".join(candidate_pid) if candidate_pid else "NA",
                                                          ",".join(mum1) if mum1 else "NA",
                                                          ",".join(dad1) if dad1 else "NA",
                                                          ",".join(mum2) if mum2 else "NA",
                                                          ",".join(dad2) if dad2 else "NA"]
                                                           )
                                comphet_ids[ped.name] += 1
    # write to disc
    comphetF=outdir+"."+chrom+".comphet_tmp"
    comphet_count=0
    with open(comphetF, 'w') as out:
        for ped_id in records_comphet.keys():
            for rdata in records_comphet[ped_id]:
                print("%s\t%s" % (ped_id, "\t".join( str(x) for x in rdata)), file=out)
                comphet_count=comphet_count+1
    del records_comphet
    tqdm.write("Finished %s with %i/%i SNV/CNV calls, %i genes and %i comphets" % (chrom, total_var_count["SNV"], total_var_count["CNV"], gene_count, comphet_count) )
    return( chrom, recF, records_genes, comphetF, ignored_SO_terms, total_var_count, filter_counts )



# .........................................................................
#                Parallel calling
# .........................................................................
try:
    queue = Manager().Queue()
    progress = Process(target=show_prog, args=(queue, total_var_count_snv + total_var_count_cnv))
    progress.start()
    # get input files
    snvF_chr_f=[]
    cnvF_chr_f=[]
    for c in chroms:
        snvF_chr_f+=[snvF_chr[c]]
        cnvF_chr_f+=[cnvF_chr[c]] if cnvF_chr is not None else None
    param = zip(
        chroms,  
        snvF_chr_f, 
        cnvF_chr_f, 
        itertools.repeat(peds), 
        itertools.repeat(gados), 
        itertools.repeat(exos), 
        itertools.repeat(gene_anno), 
        itertools.repeat(gene_anno_columns), 
        itertools.repeat(max_var_per_cat), 
        itertools.repeat(include_na_models), 
        itertools.repeat(chunksize),
        itertools.repeat(debug),
        itertools.repeat(outdir),
        itertools.repeat(queue)
        )
    with Pool(processes=threads) as pool:
        results = pool.map(process_chrom, param)
        print("receiving results...")
        progress.terminate()
        data={}
        for c in results:
            data[c[0]]= c
except:
    print("Terminating progess thread.")
    progress.terminate()
    raise       
    
     
# =================================================================================
# merge results and write result file
# =================================================================================
print("Finished, writing results.")

# =================================================================================
# Load all variant record
# =================================================================================
print("Loading output records.")
records={}
for chr in chroms:
    records[chr]={}
    chrom, recF, records_gene, comphetF, chr_ignored_SO_terms, var_count, filter_counts = data[chr]
    for ped_id in recF.keys():
        if os.stat(recF[ped_id]).st_size > 0:
            d = pd.read_csv(recF[ped_id],delimiter='\t',encoding='utf-8', header=None, dtype=str, na_values=[''], keep_default_na=False)
            for index, row in d.iterrows():
                if not ped_id in records[chr]:
                    records[chr][ped_id]=list()
                dat=[str(x) for x in row]
                records[chr][ped_id].append(dat)
        removeFile(recF[ped_id])

print("Writing variant records.")    
# for storing comphet candidates
ped2gid2vid2mumordad2pid={}
# SO terms that were ignored
ignored_SO_terms = set()
# for mapping rec_ids
rec_id_map = {}
var_id_map = {}
rec_ids = {} 
var_ids = {} 
for idx, ped in enumerate(peds):
    outF=outdir + config["dataset_name"]+ "." + ped.name + ".v2r.vars.tsv"
    with open(outF, 'w') as out:
        rec_ids[ped.name] = 0
        rec_id_map[ped.name] = {}
        var_ids[ped.name] = 0
        var_id_map[ped.name] = {}
        
        # write header
        header  = "# Samples (Affected are indicated by [*]):\n"
        for pid in ped.all_ids():
            if pid in ped.affected_ids():
                pid+=" [*]"
            header+="#    "+pid + "\n"        
        header += "# Genotype codes: hom_ref = 0, het = 1, hom_alt = 2 (which alt is untracked), uncalled = NA\n"
        header += "# Supported inheritance models are not filtered for NA genotypes or genotype quality\n"
        header += "rec_id\tvar_id\tgene\tchr\tstart\tend\tref\talt\tvar_type\taa_change\tknown_ids\tmax_pop_af\tcohort_af\treg_id\tdb_source\treg_type\tclosest_gene_dist\tconsequence\td_score\ta_score\tsup_rec\tsup_dom\tsup_dnm\tsup_comphet_mum\tsup_comphet_dad\tnum_hq_calls\thom_aff\thom_unaff\thet_aff\thet_unaff"
        for pid in ped.all_ids():
            header+="\tGT_"+pid        
        for pid in ped.all_ids():
            header+="\tGQ_"+pid        
        header += ("\t" + "\t".join(config["output_fields"]["included_info_fields"])  ) if len(config["output_fields"]["included_info_fields"])>0 else ""
        header += ("\t" + "\t".join(config["output_fields"]["included_info_factors"]) ) if len(config["output_fields"]["included_info_factors"])>0 else ""
        header += ("\t" + "\t".join(config["output_fields"]["included_info_flags"])   ) if len(config["output_fields"]["included_info_flags"])>0 else ""
        print(header, file=out)
    
        for chr in chroms:
            chrom, refF, records_gene, comphetF, chr_ignored_SO_terms, var_count, filter_counts = data[chr]
            ignored_SO_terms.update(chr_ignored_SO_terms)
            if ped.name in records[chr]:
                for rec in records[chr][ped.name]:
                    # record ids
                    rec_id_map[ped.name][rec[0]] = "r" + str(rec_ids[ped.name])
                    rec[0] = rec_id_map[ped.name][rec[0]]
                    rec_ids[ped.name] += 1
                    # cvariant ids
                    if rec[1] in var_id_map[ped.name]:
                        rec[1] = var_id_map[ped.name][rec[1]]
                    else:
                        var_id_map[ped.name][rec[1]] =  "v" + str(var_ids[ped.name])
                        rec[1] = var_id_map[ped.name][rec[1]]
                        var_ids[ped.name] += 1
                   
                    print("\t".join( str(x) for x in rec), file=out)
del records


# =================================================================================
# write genes table
# =================================================================================
print("Writing gene records.")
gene_ids = {} 
unique_genes = {}
for idx, ped in enumerate(peds):
    gene_ids[ped.name]=0
    unique_genes[ped.name]=set()
    outF=outdir + config["dataset_name"] + "." + ped.name + ".v2r.genes.tsv"
    with open(outF, 'w') as out:
        header = "# "+ped.name+"\n"
        header += "rec_id\tgene\tgado_zscore\texomiser_gene_pheno_score\tinh_model\tvariants_n\tvariants"
        if gene_anno_columns:
            header+="\t"+"\t".join(gene_anno_columns)
        print(header, file=out)
        for chr in chroms:
            chrom, records, records_gene, comphetF, chr_ignored_SO_terms, var_count, filter_counts = data[chr]
            if ped.name in records_gene:
                for rec in records_gene[ped.name]:
                    linked_vars = ",".join(rec_id_map[ped.name][str(x)] for x in rec[6])
                    rec[6] = linked_vars
                    rec[0] = "g" + str(gene_ids[ped.name])
                    print("\t".join( str(x) for x in rec), file=out)
                    gene_ids[ped.name] += 1  
                    unique_genes[ped.name].update([rec[1]])        
# =================================================================================
# write comphet table
# =================================================================================
print("Writing comphet records.")
records_comphet={}
for chr in chroms:
    records_comphet[chr]={}
    chrom, refF, records_gene, comphetF, chr_ignored_SO_terms, var_count, filter_counts = data[chr]
    if os.stat(comphetF).st_size > 0:
        d = pd.read_csv(comphetF,delimiter='\t',encoding='utf-8', header=None, dtype=str, na_values=[''], keep_default_na=False)
        for index, row in d.iterrows():
            ped_id=row[0]
            if not ped_id in records_comphet[chr]:
                records_comphet[chr][ped_id]=list()
            dat=[str(x) for x in row]
            records_comphet[chr][ped_id].append(dat[1:])
    removeFile(comphetF)
# write merged 
comphet_ids = {} 
for idx, ped in enumerate(peds):
    comphet_ids[ped.name]=0
    outF=outdir + config["dataset_name"] + "." + ped.name + ".v2r.comphet.tsv"
    with open(outF, 'w') as out:
        header = "# "+ped.name+"\n"
        header += "rec_id\tgene\tv1\tv2\tnum_aff\tcandidate\tmum1\tdad1\tmum2\tdad2"
        print(header, file=out)
        for chr in chroms:
            chrom, records, records_gene, comphetF, chr_ignored_SO_terms, var_count, filter_counts = data[chr]
            if ped.name in records_comphet[chr]:
                for rec in records_comphet[chr][ped.name]:
                    rec[0] = "c" + str(comphet_ids[ped.name])
                    rec[2] = rec_id_map[ped.name][rec[2]]
                    rec[3] = rec_id_map[ped.name][rec[3]]
                    print("\t".join( str(x) for x in rec), file=out)
                    comphet_ids[ped.name] += 1 
del records_comphet

# write index
print("Writing index.")
outF=outdir + config["dataset_name"] + ".v2r.idx.tsv"
with open(outF, 'w') as out:
    print("# V2 analysis date: "+time.strftime('%Y-%m-%d %H:%M:%S'), file=out)
    print("rec_id\tpedigree\tgene_records\tunique_genes\tvariant_records\tvariants\tcomphet_records\tn_all_samples\tn_affected_samples\tall_samples\taffected_samples\tgene_file\tvariant_file\tcomphet_file\tpedigree_file\tgado_file\texomiser_files", file=out)
    rec_id = 0
    for idx, ped in enumerate(peds):
        genf = os.path.abspath(outdir + config["dataset_name"] + "." + ped.name + ".v2r.genes.tsv")
        bgzip(genf, index=True, delinFile=True, maxthreads=threads)
        genf+=".gz"
        
        varf = os.path.abspath(outdir + config["dataset_name"] + "." + ped.name + ".v2r.vars.tsv")
        bgzip(varf, index=True, delinFile=True, maxthreads=threads)
        varf+=".gz"
        
        comf = os.path.abspath(outdir + config["dataset_name"] + "." + ped.name + ".v2r.comphet.tsv")
        bgzip(comf, index=True, delinFile=True, maxthreads=threads)
        comf+=".gz"
        
        pedf = os.path.abspath(config["input_data"]["ped_pattern"].format(id=ped.name))
        gadof = os.path.abspath(config["input_data"]["gado_pattern"].format(id=ped.name))
        exofs=[]
        for im in EXOMISER_IM:
            exof = config["input_data"]["exomiser_pattern"].format(IM=im, id=ped.name)
            if files_exist(exof):
                exofs+=[exof]
        rec=[rec_id,
            ped.name, 
            gene_ids[ped.name], 
            len(unique_genes[ped.name]),
            rec_ids[ped.name],        
            var_ids[ped.name],         
            comphet_ids[ped.name],             
            len(ped.all_ids()), len(ped.affected_ids()), 
            ",".join(ped.all_ids()),",".join(ped.affected_ids()),
            genf, varf, comf, pedf, gadof, ",".join(exofs)]
        print("\t".join( str(x) for x in rec), file=out)
        rec_id += 1
bgzip(outF, index=True, delinFile=True, maxthreads=threads)
    

# write merged debug log
dfiles = [outdir+"."+chr+".log_tmp" for chr in chroms]
if debug:
    with open(outdir+"var2reg.debug.log", 'w') as outfile:
        for fname in dfiles:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
removeFile(dfiles)
                
with open(outdir+"var2reg.stats.tsv", 'w') as log:
    # write stats
    print("#=============================================", file=log)
    print("# Filtered variants statistics", file=log)
    print("# AF = maxPopAF too high", file=log)
    print("# NC = no call in any affected sample", file=log)
    print("# LS = low score", file=log)
    print("# NI = no supported inheritance model (only if include_na_models is False)", file=log)
    print("# ID = Incomplete data (CNVs only)", file=log)
    print("# UT = unspecific variant type (CNVs only)", file=log)
    print("# NP = no_pass (these variants are not included in the debug log!)", file=log)
    print("=============================================", file=log)
    all_var_count={}
    all_filter_sum={}
    all_filter_counts={}
    for vcf_type in ["SNV", "CNV"]:
        all_filter_counts[vcf_type]={}
        all_var_count[vcf_type]=0
        all_filter_sum[vcf_type]=0
        for tag in filter_tags:
            all_filter_counts[vcf_type][tag]=0
    print("var_type\tChr\tVariants\tfrac_filtered\t%s\t%s" % ("\t".join(filter_tags), "\t".join("frac_"+x for x in filter_tags) ), file=log )
    for vcf_type in ["SNV", "CNV"]:
        for chr in chroms:
            chrom, refF, records_gene, comphetF, chr_ignored_SO_terms, var_count, filter_counts = data[chr]
            fsum = sum(filter_counts[vcf_type].values())
            all_var_count[vcf_type] += var_count[vcf_type]
            all_filter_sum[vcf_type] += fsum
            for tag in filter_tags:
                all_filter_counts[vcf_type][tag]+=filter_counts[vcf_type][tag]
            print("%s\t%s\t%i\t%f\t%s\t%s" % (vcf_type, 
                                              chrom, 
                                              var_count[vcf_type], 
                                              float(fsum)/var_count[vcf_type] if var_count[vcf_type]>0 else 0, 
                                              "\t".join(str(filter_counts[vcf_type][x]) for x in filter_tags), 
                                              "\t".join(str(filter_counts[vcf_type][x]/float(fsum) if fsum > 0 else 0) for x in filter_tags) ), file=log)
        print("%s\t%s\t%i\t%f\t%s\t%s" % (vcf_type, "ALL",
                                  all_var_count[vcf_type], 
                                  float(all_filter_sum[vcf_type])/all_var_count[vcf_type] if all_var_count[vcf_type]>0 else 0, 
                                  "\t".join(str(all_filter_counts[vcf_type][x]) for x in filter_tags), 
                                  "\t".join(str(all_filter_counts[vcf_type][x]/float(all_filter_sum[vcf_type]) if all_filter_sum[vcf_type] > 0 else 0) for x in filter_tags) ), file=log)
    total_var_count=0
    total_filter_sum=0
    total_filter_counts={}
    for tag in filter_tags:
        total_filter_counts[tag]=0
    for vcf_type in ["SNV", "CNV"]:
        total_var_count+=all_var_count[vcf_type]
        total_filter_sum+=all_filter_sum[vcf_type]
        for tag in filter_tags:
            total_filter_counts[tag]+=all_filter_counts[vcf_type][tag]
    print("%s\t%s\t%i\t%f\t%s\t%s" % ("ALL", "ALL",
                              total_var_count, 
                              float(total_filter_sum)/total_var_count if total_var_count>0 else 0, 
                              "\t".join(str(total_filter_counts[x]) for x in filter_tags), 
                              "\t".join(str(total_filter_counts[x]/float(total_filter_sum) if total_filter_sum > 0 else 0) for x in filter_tags) ), file=log)

# =================================================================================
# clean up
# =================================================================================
removeFile(snvF_chr_f)
removeFile([x+".tbi" for x in snvF_chr_f])
removeFile(cnvF_chr_f)
removeFile([x+".tbi" for x in cnvF_chr_f])

print("\nNOTE: the following SO terms were ignored: %s" % (ignored_SO_terms))
print("Finished in " + str(datetime.timedelta(seconds=time.time()-startTime)))
print("Done.")
