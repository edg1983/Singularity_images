
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

from gffutils.iterators import DataIterator
from numpy.distutils.fcompiler import none
import progressbar
from vcfpy import *
import vcfpy

import pandas as pd
import pyranges as pr
import warnings

# Necessary for including python modules from a parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
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
"""


# @see https://stackoverflow.com/questions/2187269/print-only-the-message-on-warnings
def custom_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return str(msg) + '\n'
warnings.formatwarning = custom_formatwarning

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

# # FIXME also report how score was calculated
# def calc_nc_dscore(rec, conf):   
#     scores=[]
#     score_cadd = rec.INFO.get("CADD_PhredScore", None) # range: 0-35
#     if score_cadd is not None:
#         scores+=[toFloat(score_cadd) / 35.0]
#     score_dann = rec.INFO.get("DANN_DANN", None) # range 0-1
#     if score_dann is not None:
#         scores+=[toFloat(score_dann)]
#     score_lin = rec.INFO.get("LinSight", None) # range 0-1
#     if score_lin is not None:
#         scores+=[toFloat(score_lin)]
#     score_remm = rec.INFO.get("ReMM_score", None) # range 0-1
#     if score_remm is not None:
#         scores+=[toFloat(score_remm)]
#     sum=0
#     for s in scores:
#         sum+=s
#     d_score = float(sum)/float(len(scores)) if len(scores)>0 else float(0)
#     return (d_score)
# 
# 
# # maximum of NC score and preconfigured weight (if any)
# def calc_intron_dscore(rec, d_score_conf, conf):
#     w_nc = calc_nc_dscore(rec, conf)
#     w_splice = float(rec.INFO.get("SpliceAI_DS_MAX", 0)) # range: 0-1  
#     return getMax([w_nc, w_splice, d_score_conf])
# 
# # maximum of NC score and preconfigured weight (if any)
# def calc_utr_dscore(rec, d_score_conf, conf):
#     w_nc = calc_nc_dscore(rec, conf)
#     return getMax([w_nc, d_score_conf])
# 
# # max of preconfigured weight (if any) and splice AI score or NC score 
# def calc_splicing_dscore(rec, d_score_conf, conf):
#     w_splice = float(rec.INFO.get("SpliceAI_DS_MAX", 0)) # range: 0-1
#     return getMax([w_splice, d_score_conf, conf])
#                  
# def calc_max_AF(rec, info_fields):
#     dat=[]
#     for f in info_fields:
#         dat.append(rec.INFO.get(f, None))
#     return getMax(dat)   


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
    def calc_inh_support(self, rec, min_gq):
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
                        if gq_m is not None and gq_m > min_gq and gq_d is not None and gq_d > min_gq:
                            if gt_m == 0 and gt_d == 0:
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


# GT : genotype, encoded as allele values separated by either of / or |. 
# The allele values are 0 for the reference allele (what is in the REF field), 1 for the first allele listed in 
# ALT, 2 for the second allele list in ALT and so on. For diploid calls examples could be 0/1, 1 | 0, or 1/2, etc. 
# For haploid calls, e.g. on Y, male nonpseudoautosomal X, or mitochondrion, only one allele value should be 
# given; a triploid call might look like 0/0/1. If a call cannot be made for a sample at a given locus, 
# ‘.’ should be specified for each missing allele 5 in the GT field (for example ‘./.’ for a diploid genotype 
# and ‘.’ for haploid genotype). The meanings of the separators are as follows (see the PS field below for 
# more details on incorporating phasing information into the genotypes): / : genotype unphased | : genotype phased
class TYPE(Enum):
    HOM = 1
    HET = 2 
    HOMREF = 3
    NOCALL = 4   


class Genotype():

    def __init__(self, str):
        super().__init__()
        self.str = str
        self.all = str.split("|") if "|" in str else str.split("/")

    def __repr__(self):
        return self.__str__()  

    def __str__(self):
        return ('[GT: %s, %s]' % (self.type(), ", ".join(self.all)))

    def is_het(self):
        return self.type() is TYPE.HET

    def is_hom(self):
        return self.type() is TYPE.HOM

    def is_homref(self):
        return self.type() is TYPE.HOMREF

    def is_nocall(self):
        return self.type() is TYPE.NOCALL

    def type(self, allele=1):
        dots = 0
        hits = 0
        for a in self.all:
            if a == ".":
                dots = dots + 1
            if a == str(allele):
                hits = hits + 1
        if dots == len(self.all):  # all dots: this is a no_call
            return TYPE.NOCALL  
        if hits == len(self.all):  # all hits: this is a HOM call
            return TYPE.HOM   
        if len(self.all) == 2 and hits == 1:  # one hit: this is a REF call (including e.g., "./1")
            return TYPE.HET   
        # treat all other situations as nocalls
        return TYPE.NOCALL  
    
    
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
    
    print(a)
    print(b)
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
parser.add_argument("-ids","--ids", type=str, required=True, nargs='+', dest="ids", metavar="ids", help="Sample ids [must match ped and gado filenames]")
parser.add_argument("-ped","--ped", type=str, required=True, dest="pedFtemplate", metavar="pedFtemplate", help="Input PED file template. possible tokens: {id}")
parser.add_argument("-gad","--gad", type=str, required=True, dest="gadoFtemplate", metavar="gadoFtemplate", help="Input GADO file template. possible tokens: {id}")
parser.add_argument("-exo","--exo", type=str, required=True, dest="exoFtemplate", metavar="exoFtemplate", help="Input Exomiser file pattern (e.g., /mydir/{id}.{IM}.genes.tsv where {id} will be replaced by the respective id and {IM} sequentially by the following inheritance model ids: (" + ",".join(EXOMISER_IM) +")")
parser.add_argument("-reg","--reg", type=existing_file, required=True, dest="regF", metavar="regF", help="Input regulatory regions")
parser.add_argument("-snv","--snv", type=existing_file, required=True, dest="snvF", metavar="snvF", help="Input SNV VCF or pre-filtered V2 TSV file")
parser.add_argument("-cnv","--cnv", type=existing_file, required=False, dest="cnvF", metavar="cnvF", help="Input CNV VCF file")
parser.add_argument("-gff","--gff", type=existing_file, required=True, dest="gffF", metavar="gffF", help="Input GFF3 annotation file (used to filter intronic variants)")
parser.add_argument("-gen","--gen", type=existing_file, required=False, dest="genF", metavar="genF", help="Input gene annotation file")
parser.add_argument("-pre","--pre", type=str, default="V2", required=False, dest="prefix", metavar="prefix", help="File name prefix")
parser.add_argument("-ali","--alias", type=existing_file, required=False, dest="aliasF", metavar="aliasF", help="Optional gene alias file")
parser.add_argument("-o","--out", type=str, required=True, dest="outdir", metavar="outdir", help="output folder")
parser.add_argument("-conf","--config", type=existing_file, required=True, dest="confF", help="JSON config file")
parser.add_argument("--write_vcf", required=False, action="store_true", default=False, dest="write_vcf", help="If set, VCF files will be written (useful for debugging)")
parser.add_argument("--include_na_models", required=False, action="store_true", default=False, dest="include_na_models", help="If set, output will also contain variants that support no inheritance model (useful for debugging)")
parser.add_argument("--include_nocalls", required=False, action="store_true", default=False, dest="include_nocalls", help="If set, output VCFs will also contain variants w/o calls in affected sample (useful for debugging)")
parser.add_argument("--include_info_fields", required=False, action="store_true", default=False, dest="include_info_fields", help="If set, output VCFs will also include input INFO fields (useful for debugging)")
parser.add_argument("--max_var_per_cat", type=int, required=False, dest="max_var_per_cat", help="debug only")

args = parser.parse_args() 


print(logo)
startTime = time.time()

# ignore warnings (mainly from vcfpy)
if not sys.warnoptions:
    warnings.filterwarnings(action="ignore", module="vcfpy")
#============================================================================
chunksize = 1000000


# ensure dirs
if not os.path.exists(args.outdir):
        print("Creating dir " + args.outdir)
        os.makedirs(args.outdir)
log = os.path.join(args.outdir, pipename.replace(" ", "-") + ".log")
if not args.outdir.endswith("/"):
    args.outdir += "/"
    
config = json.load(open(args.confF), object_pairs_hook=OrderedDict)
# add commandline params to config
config["CMD"]={}
for arg in vars(args):
    config["CMD"][arg] = getattr(args, arg) 
# write effective config
outF=args.outdir + args.prefix + ".var2reg.effective_conf.json"
with open(outF, 'w') as outs:
    print(json.dumps(config, indent=4, sort_keys=False), file=outs)

chroms = None
if "chrom" in config:
    chroms = config["chrom"]
    print("Handling chromosomes %s" % (", ".join(chroms)))

# read gene name alias table
alias2gid={}
if args.aliasF:
    d = pd.read_csv(args.aliasF,delimiter='\t',encoding='utf-8')
    alias2gid=pd.Series(d.symbol.values,index=d.alias).to_dict()
print("Loaded %i gene symbol aliases" % (len(alias2gid)))


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
if args.genF:
    nline = sum(1 for i in gzip.open(args.genF, 'rb'))
    print("Accessing file with %i additional gene annotations" % (nline))
    with progressbar.ProgressBar(max_value=nline) as bar:  
        for d in pd.read_csv(args.genF,
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
nline = sum(1 for i in gzip.open(args.regF, 'rb'))
print("Accessing file with %i regulatory regions" % (nline))
with progressbar.ProgressBar(max_value=nline) as bar:  
    names = ['Chromosome', 'Start', 'End', 'id', 'std_type', 'db_source', 'closestGene_symbol', 'closestGene_dist', 'affected_genes', 'N_methods'] 
    dtype = {'Chromosome': str, 'Start': int, 'End': int, 'id': str, 'std_type': str, 'db_source': str, 'closestGene_symbol': str, 'closestGene_dist': int, 'affected_genes': str, 'N_methods': str}
    reg_regions=[]
    id2reg_regions={}
    for d in pd.read_csv(args.regF,
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
gff = pr.read_gff3(args.gffF)
exon_int = gff[gff.Feature == 'exon']
utr3_int = gff[gff.Feature == 'three_prime_UTR']
utr5_int = gff[gff.Feature == 'five_prime_UTR']
utr_int = pr.concat([utr3_int, utr5_int])
gene_int = gff[gff.Feature == 'gene']
print("Loaded %i exon, %i utr and %i gene annotations from %s" % (len(exon_int), len(utr_int), len(gene_int), args.gffF))

# read samples we have data for from VCF header
reader = vcfpy.Reader.from_path(args.snvF)
reader2 = vcfpy.Reader.from_path(args.cnvF) if args.cnvF else None


# =================================================================================
# read input pedigrees
# =================================================================================
print("Reading pedigrees and gene rank data")
peds=[]
gados=[]
exos=[]
for id in args.ids:
    pedf = args.pedFtemplate.format(id=id)
    gadof = args.gadoFtemplate.format(id=id)
    anyexo = False
    for im in EXOMISER_IM:
        exof = args.exoFtemplate.format(IM=im, id=id)
        anyexo = anyexo or files_exist(exof)
        #print(args.exoDir + "/" + id + "_WGS_ALL_" + im + ".genes.tsv")
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
                exof = args.exoFtemplate.format(IM=im, id=id)
                if files_exist([exof]):
                    d = pd.read_csv(exof, delimiter='\t', encoding='utf-8')
                    for index, row in d.iterrows():
                        exo2rnk[row['#GENE_SYMBOL']]=row['EXOMISER_GENE_PHENO_SCORE']
            exos+=[exo2rnk]       
            print("\t%s: %i gado_rnk, %i exo_ps" % (ped.name, len(gid2rnk), len(exo2rnk)))
        else:
            print("\tDid not find any data for the following samples in the snv VCF: %s. Ignoring ped %s." % (",".join(ped.all_ids()), id))
    else:
        print("\tCannot load all required data for id %s. Some files (ped, gado or exomiser results) do not exist. Ignoring (%s, %s, %s)" % (id, pedf, gadof, args.exoFtemplate) )
    #print("\tsamples: %s" % (ped.all_ids()))
    #print("\taffected_samples: %s" % (ped.affected_ids()))
print("\tDone. Read %i pedigrees and %i gado rank lists" % (len(peds), len(gados)))
if len(peds) == 0:
    sys.exit("Check input configuration / data...")

# =================================================================================
# prepare
# =================================================================================
outs = {}
out_ids = {} # counter of output records (there maybe multiple record per variant)
out_vcf = {}
for vcf_type in ["SNV", "CNV"]:
    out_vcf[vcf_type] = {}
for idx, ped in enumerate(peds):
    outF=args.outdir + args.prefix + "." + ped.name + ".var2reg.vars.tsv"
    outs[ped.name] = open(outF, 'w')
    out_ids[ped.name] = 0 
    # vcf output
    outVcfSnvF=args.outdir + args.prefix + "." + ped.name + ".var2reg.vars.snv.vcf"
    header_snv = copy.deepcopy(reader.header)
    header_snv.samples = vcfpy.SamplesInfos(sample_names= [x for x in ped.all_ids() if x in reader.header.samples.names]  )
    if args.write_vcf:
        out_vcf["SNV"][ped.name] = vcfpy.Writer.from_path(outVcfSnvF, header_snv)
    if args.cnvF:
        outVcfCnvF=args.outdir + args.prefix + "." + ped.name + ".var2reg.vars.cnv.vcf"
        header_cnv = copy.deepcopy(reader2.header)
        header_cnv.samples = vcfpy.SamplesInfos(sample_names=[x for x in ped.all_ids() if x in reader2.header.samples.names])
        if args.write_vcf:
            out_vcf["CNV"][ped.name] = vcfpy.Writer.from_path(outVcfCnvF, header_cnv)
    # write header
    header  = "# Samples (Affected are indicated by [*]):\n"
    for pid in ped.all_ids():
        if pid in ped.affected_ids():
            pid+=" [*]"
        header+="#    "+pid + "\n"        
    header += "# Genotype codes: hom_ref = 0, het = 1, hom_alt = 2 (which alt is untracked), uncalled = NA\n"
    header += "# Supported inheritance models are not filtered for NA genotypes or genotype quality\n"
    header += "ID\tVar2reg_id\tGene\tChr\tPos\tEnd\tRef\tAlt\tVarType\tknownIds\tMaxPopAF\tcohortAF\tRegion_id\tDb_source\tReg_type\tClosest_gene_dist\tConsequence\td_score\ta_score\tsup_rec\tsup_dom\tsup_dnm\tsup_comphet_mum\tsup_comphet_dad\tnum_hq_calls\thom_aff\thom_unaff\thet_aff\thet_unaff"
    for pid in ped.all_ids():
        header+="\tGT_"+pid        
    for pid in ped.all_ids():
        header+="\tGQ_"+pid        
    header += ("\t" + "\t".join(config["included_info_fields"])  ) if len(config["included_info_fields"])>0 else ""
    header += ("\t" + "\t".join(config["included_info_factors"]) ) if len(config["included_info_factors"])>0 else ""
    header += ("\t" + "\t".join(config["included_info_flags"])   ) if len(config["included_info_flags"])>0 else ""
    print(header, file=outs[ped.name])


# =================================================================================
# iterate variants and write var2reg tables
# =================================================================================
# for storing vars
ped2gid2mod2vid={}
# for storing comphet candidates
ped2gid2vid2mumordad2pid={}
# SO terms that were ignored
ignored_SO_terms = set()
# for keeping track of current variant
wrote2vcf={} 
# variant id
var_id=0
for vcf_type in ["SNV", "CNV"]:
    print("Analyzing %s records" % (vcf_type))
    # for warning about slow pyranges queries
    did_warn_pyranges=False

    # load data
    infile = args.snvF if vcf_type == "SNV" else args.cnvF
    if not infile:
        print("No VCF file configured for type %s." % infile )
        continue # no CNV file configured
    nline = sum(1 for i in gzip.open(infile, 'rb'))
    print("Parsing ~%i %s from %s" % (nline, vcf_type, infile))
    # read VCF header via vcfpy (# full file: 13Mio SNVs)
    header = reader.header if vcf_type=="SNV" else reader2.header
    parser = vcfpy.parser.RecordParser(header, header.samples)
    if vcf_type=="CNV":
            # fix issue with missing CN header in CNV vcf
            header.add_format_line(vcfpy.OrderedDict([('ID', 'CN'), ('Description', 'Copynumber'), ('Number', '1'), ('Type', 'Float')]))
    names = ['CHROM', 'POS', 'ID', 'REF', 'ALT',
               'QUAL', 'FILTER', 'INFO', 'FORMAT'] 
    dtype = {'CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
       'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT':str}
    for s in header.samples.names:
        names = names + [s]
        dtype[s] = str
    with progressbar.ProgressBar(max_value=nline) as bar:  

        var_count = 0
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
            # debug
            if args.max_var_per_cat and var_calls >= args.max_var_per_cat:
                break
        
            # pass only
            d = d.query("FILTER=='PASS' | FILTER=='.'")
            if chroms:
                d = d[d['CHROM'].isin(chroms)]
             
            for index, row in d.iterrows():
                var_count += 1
                if var_count % 10 == 0:
                    bar.update(var_count)
                    # debug
                if args.max_var_per_cat and var_calls >= args.max_var_per_cat:
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
                    if args.write_vcf:
                        out_vcf[vcf_type][ped.name].write_record(get_output_rec(rec, args.include_info_fields,'INCOMPLETE_DATA',None) )
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
                        if args.write_vcf:
                            out_vcf[vcf_type][ped.name].write_record(get_output_rec(rec, args.include_info_fields,'UNSPECIFIC_TYPE',None) )
                        continue
                    # this will map, e.g., DEL=>DEL, DEL:ME=>DEL:ME, DEL:ME:something=>DEL:ME
                    vartype = ":".join(vartype.split(":")[0:2])
                    
                # filter by max pop AF
                #maxAF = calc_max_AF(rec, ["global_A1000G_AF", "global_ExAC_AF", "global_UK10K_AF", "global_gnomAD_AF"] if vcf_type =="SNV" else ["GD_POPMAX_AF", "1000g_max_AF"])
                maxAF = calc_score(rec, config, "max_pop_af_snv", 0) if vcf_type =="SNV" else calc_score(rec, config,"max_pop_af_cnv",0)
                if maxAF is not None and maxAF > config["max_pop_af"]:
                    info_add={}
                    info_add["maxAF"]=str(maxAF)
                    if args.write_vcf:
                        out_vcf[vcf_type][ped.name].write_record(get_output_rec(rec, args.include_info_fields,'TOO_COMMON',None, info_add=info_add) )
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
                    if gt is not None and maxgq > config["min_gq"]:
                        cohort_calls += gt
                        all_calls += 2
                cohortAF = 0 if all_calls == 0 else float(cohort_calls) / float(all_calls) # allele frequency in cohort
                
                # iterate over all pedigrees
                for idx, ped in enumerate(peds):
                    out = outs[ped.name]
                    wrote2vcf[ped.name]=False
                
                    # calc support for inh models
                    min_gq = config["min_gq"] if vcf_type == "SNV" else -1 # GQ is not a reliable measure for CNVs
                    sup_rec, sup_dom, sup_dnm, num_hq_calls, hom_aff, hom_unaff, het_aff, het_unaff, gt_str, gq_str, sup_comhet_mum, sup_comhet_dad = ped.calc_inh_support(rec, min_gq)
                    max_sup = getMax([sup_rec, sup_dom, sup_dnm])
                    sup_comhet = len(sup_comhet_mum) + len(sup_comhet_dad)
                    num_affected_calls = hom_aff + het_aff
                    
                    # no call with GQ>thresh in affected samples. ignore.
                    if num_affected_calls == 0: 
                        if args.include_nocalls and args.write_vcf:
                            out_vcf[vcf_type][ped.name].write_record(get_output_rec(rec, args.include_info_fields, 'NO_CALL', None))
                        continue
                    # no support for any inh model or comphet. ignore.
                    if not args.include_na_models and sup_comhet == 0 and ( max_sup is None or max_sup <= 0):
                        if args.write_vcf:
                            info_add={}
                            info_add["max_sup"]=str(max_sup)
                            out_vcf[vcf_type][ped.name].write_record(get_output_rec(rec, args.include_info_fields,'NO_INH_MODEL',None, info_add=info_add) )
                        continue
        
                    # collect INFO values   
                    infos=[]
                    for f in config["included_info_fields"]:
                        infos.append(fmt(getMax(str(rec.INFO.get(f, "NA")).split(","))))
                    for f in config["included_info_factors"]:
                        infos.append(str(rec.INFO.get(f, "NA")))
                    for f in config["included_info_flags"]:
                        infos.append("1" if f in rec.INFO else "0")
                        
                    #--------------------------------
                    # Find overlap with genic regions
                    #--------------------------------
                    # get term with max weight per gene
                    gene2term={}
                    if vcf_type == "SNV":
                        if "ANN" in rec.INFO:
                            for a in rec.INFO["ANN"]:
                                for aa in rec.ALT:
                                    dat = a.split("|")
                                    gid=dat[3]
                                    gid = alias2gid[gid] if gid in alias2gid else gid
                                    #print(dat[0]+"/"+aa.serialize() + "-" + dat[1])
                                    if dat[0] == aa.serialize():
                                        for imp in dat[1].split("&"):
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
                        
                        if vcf_type == "CNV" or d_score >= config["min_dscore"]:
                            handled_gids += gid
                            print("%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%i\t%i\t%i\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%s\t%s\t%s" %
                                (out_ids[ped.name], var_id, gid, rec.CHROM, pos, end, rec.REF, 
                                altstr, vartype, ",".join(rec.ID) if len(rec.ID)>0 else "NA", fmt(maxAF), fmt(cohortAF), 
                                "NA", "anno", term.type, "NA", term.term, fmt(d_score), fmt(a_score), 
                                sup_rec, sup_dom, sup_dnm, 
                                ",".join(sup_comhet_mum) if sup_comhet_mum else "NA", ",".join(sup_comhet_dad) if sup_comhet_dad else "NA",  
                                num_hq_calls, hom_aff, hom_unaff, het_aff, het_unaff,
                                gt_str, gq_str, "\t".join(infos)), file=out)
        
                            # add gene<>inh_model<>var links
                            if sup_rec > 0:
                                ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "recessive", out_ids[ped.name])
                            if sup_dom > 0:
                                ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "dominant", out_ids[ped.name])
                            if sup_dnm > 0:
                                ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "denovo", out_ids[ped.name])
                            if args.include_na_models:
                                ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "NA", out_ids[ped.name])
        
                            # comphet. sup_comhet_mum, sup_comhet_dad contain the pids of affected samples that inherit a het from mum/dad
                            for pid in sup_comhet_mum:
                                ped2gid2vid2mumordad2pid = add_link2(ped2gid2vid2mumordad2pid, ped.name, gid, out_ids[ped.name], "mum", pid)
                            for pid in sup_comhet_dad:
                                ped2gid2vid2mumordad2pid = add_link2(ped2gid2vid2mumordad2pid, ped.name, gid, out_ids[ped.name], "dad", pid)
        
                            # inc var id and write VCF
                            if not wrote2vcf[ped.name] and args.write_vcf:
                                out_vcf[vcf_type][ped.name].write_record(get_output_rec(rec, args.include_info_fields,None,str(var_id)) )
                            wrote2vcf[ped.name]=True
                            out_ids[ped.name] += 1
                            var_calls += 1
                        
                    #--------------------------------
                    # Find overlap with regulatory regions. Filter for dscore only for SNVs
                    #--------------------------------    
                    if vcf_type == "SNV":
                        d_score=calc_score(rec, config, "noncoding", 0)
                    else:
                        d_score=None
                    if vcf_type == "CNV" or d_score >= config["min_dscore"]:
                        gid2reg={}
                        gid2assoctype={}
                        if vcf_type == "SNV" and "snv_reg_id_info_field" in config:
                            # try to get reg_id values from INFO fields
                            if rec.INFO.get(config["snv_reg_id_info_field"], None):
                                for rid in rec.INFO.get(config["snv_reg_id_info_field"]).split(","):
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
                        elif vcf_type == "CNV" and "cnv_reg_id_info_field" in config:
                            # try to get reg_id values from INFO fields
                            if rec.INFO.get(config["cnv_reg_id_info_field"], None):
                                for rid in rec.INFO.get(config["snv_reg_id_info_field"]).split(","):
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
                                print("%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%i\t%i\t%i\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%s\t%s\t%s" %
                                    (out_ids[ped.name], var_id, gid, rec.CHROM, pos, end, rec.REF, 
                                     altstr, vartype, ",".join(rec.ID) if len(rec.ID)>0 else "NA", fmt(maxAF), fmt(cohortAF), 
                                     region['id'], region['db_source'], gid2assoctype[gid], region['closestGene_dist'], region['std_type']+"_variant", 
                                     fmt(d_score), fmt(a_score), 
                                     sup_rec, sup_dom, sup_dnm,
                                     ",".join(sup_comhet_mum) if sup_comhet_mum else "NA", ",".join(sup_comhet_dad) if sup_comhet_dad else "NA",
                                     num_hq_calls, hom_aff, hom_unaff, het_aff, het_unaff,
                                     gt_str, gq_str, "\t".join(infos) ), file=out)
                                # add gene<>inh_model<>var links
                                if sup_rec > 0:
                                    ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "recessive", out_ids[ped.name])
                                if sup_dom > 0:
                                    ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "dominant", out_ids[ped.name])
                                if sup_dnm > 0:
                                    ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "denovo", out_ids[ped.name])
                                if args.include_na_models:
                                    ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "NA", out_ids[ped.name])
        
                                # comphet. sup_comhet_mum, sup_comhet_dad contain the pids of affected samples that inherit a het from mum/dad
                                for pid in sup_comhet_mum:
                                    ped2gid2vid2mumordad2pid = add_link2(ped2gid2vid2mumordad2pid, ped.name, gid, out_ids[ped.name], "mum", pid)
                                for pid in sup_comhet_dad:
                                    ped2gid2vid2mumordad2pid = add_link2(ped2gid2vid2mumordad2pid, ped.name, gid, out_ids[ped.name], "dad", pid)
        
                                # inc var id and write VCF
                                if not wrote2vcf[ped.name] and args.write_vcf:
                                    out_vcf[vcf_type][ped.name].write_record(get_output_rec(rec, args.include_info_fields,None,str(var_id)) )
                                out_ids[ped.name] += 1
                                var_calls +=1
                                wrote2vcf[ped.name]=True
                    # write filtered
                    if not wrote2vcf[ped.name]:
                        if args.write_vcf:
                            out_vcf[vcf_type][ped.name].write_record(get_output_rec(rec, args.include_info_fields,'LOW_SCORE',None) )
                    else:
                        var_id+=1        
                    
        bar.update(nline)
    print("Wrote %i %s calls" % (var_id, vcf_type))


# close output streams
for ped in peds:
    outs[ped.name].close()
    for vcf_type in ["SNV", "CNV"]:
        if ped.name in out_vcf[vcf_type] and args.write_vcf:
            out_vcf[vcf_type][ped.name].close()
    

# write gene tables
gene_ids = {} 
for idx, ped in enumerate(peds):
    outF=args.outdir + args.prefix + "." + ped.name + ".var2reg.genes.tsv"
    with open(outF, 'w') as outs:
        header = "# "+ped.name+"\n"
        header += "Var2reg_id\tGene\tGado_zscore\tExomiser_GenePhenoScore\tInh_model\tVariants_n\tVariants"
        if gene_anno_columns:
            header+="\t"+"\t".join(gene_anno_columns)
        print(header, file=outs)
        gene_ids[ped.name]=0
        if ped.name in ped2gid2mod2vid:
            for gid in ped2gid2mod2vid[ped.name]:
                for ms in ped2gid2mod2vid[ped.name][gid]:
                    gado_z = str(gados[idx][gid]) if gid in gados[idx] else "NA"
                    exo_ps = str(exos[idx][gid]) if gid in exos[idx] else "NA"
                    annostr=""
                    if gene_anno_columns:
                        anno=gene_anno[gid] if gid in gene_anno else ["NA"] * len(gene_anno_columns)
                        annostr="\t" + "\t".join(str(x) for x in anno)
                    print("%i\t%s\t%s\t%s\t%s\t%i\t%s%s" % (gene_ids[ped.name], gid, gado_z, exo_ps, ms, len(ped2gid2mod2vid[ped.name][gid][ms]), ",".join(str(x) for x in ped2gid2mod2vid[ped.name][gid][ms]), annostr),file=outs)
                gene_ids[ped.name] += 1
                

# write comphet table
comphet_ids = {} 
for idx, ped in enumerate(peds):
    outF=args.outdir + args.prefix + "." + ped.name + ".var2reg.comphet.tsv"
    with open(outF, 'w') as outs:
        header = "# "+ped.name+"\n"
        header += "Var2reg_id\tGene\tV1\tV2\tnum_aff\tcandidate\tmum1\tdad1\tmum2\tdad2"
        print(header, file=outs)
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
                                print("%i\t%s\t%i\t%i\t%i\t%s\t%s\t%s\t%s\t%s" % (comphet_ids[ped.name], 
                                                          gid, 
                                                          v1,
                                                          v2,
                                                          len(candidate_pid), 
                                                          ",".join(candidate_pid) if candidate_pid else "NA",
                                                          ",".join(mum1) if mum1 else "NA",
                                                          ",".join(dad1) if dad1 else "NA",
                                                          ",".join(mum2) if mum2 else "NA",
                                                          ",".join(dad2) if dad2 else "NA",
                                                           ),file=outs)
                                comphet_ids[ped.name] += 1
                

# write index
outF=args.outdir + args.prefix + ".var2reg.idx.tsv"
with open(outF, 'w') as outs:
    print("# V2 analysis date: "+time.strftime('%Y-%m-%d %H:%M:%S'), file=outs)
    print("ID\tUnique_genes\tVariant_records\tComphet_records\tn_AllSamples\tn_AffectedSamples\tAllSamples\tAffectedSamples\tGeneF\tVariantF\tComphetF\tPedigreeF\tGadoF\tExoF", file=outs)
    for idx, ped in enumerate(peds):
        genf = os.path.abspath(args.outdir + args.prefix + "." + ped.name + ".var2reg.genes.tsv")
        varf = os.path.abspath(args.outdir + args.prefix + "." + ped.name + ".var2reg.vars.tsv")
        comf = os.path.abspath(args.outdir + args.prefix + "." + ped.name + ".var2reg.comphet.tsv")
        pedf = os.path.abspath(args.pedFtemplate.format(id=ped.name))
        gadof = os.path.abspath(args.gadoFtemplate.format(id=ped.name))
        exofs=[]
        for im in EXOMISER_IM:
            exof = args.exoFtemplate.format(IM=im, id=ped.name)
            if files_exist(exof):
                exofs+=[exof]
        
        print("%s\t%i\t%i\t%i\t%i\t%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
            ped.name, 
            gene_ids[ped.name], 
            out_ids[ped.name],         
            comphet_ids[ped.name],             
            len(ped.all_ids()), len(ped.affected_ids()), 
            ",".join(ped.all_ids()),",".join(ped.affected_ids()),
            genf, varf, comf, pedf, gadof, ",".join(exofs))
            , file=outs)
        
    
print("\nNOTE: the following SO terms were ignored: %s" % (ignored_SO_terms))
print("Finished in " + str(datetime.timedelta(seconds=time.time()-startTime)))
print("Done.")
