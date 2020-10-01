#!/usr/bin/env python3
'''
@author: niko
'''

from argparse import ArgumentParser, RawDescriptionHelpFormatter, \
    ArgumentTypeError
import collections
import csv, datetime, time, logging, sys, os, json
from enum import Enum
import gzip
import re
import math
from subprocess import *

from gffutils.iterators import DataIterator
from numpy.distutils.fcompiler import none
import progressbar
from vcfpy import *
import vcfpy

import pandas as pd
import pyranges as pr

# Necessary for including python modules from a parent directory
sys.path.append("/opt/V2_core")
from core.utils import *

#============================================================================

pipename = "V2 comphet calling pipeline"

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

# see https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
HIGH_IMPACT = ["transcript_ablation",  # A feature ablation whereby the deleted region includes a transcript feature
               "splice_acceptor_variant",  # A splice variant that changes the 2 base region at the 3' end of an intron
               "splice_donor_variant",  # A splice variant that changes the 2 base region at the 5' end of an intron
               "stop_gained",  # A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript
               "frameshift_variant",  # A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three
               "stop_lost",  # A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript
               "transcript_amplification",  # A codon variant that changes at least one base of the canonical start codon
               "start_lost",  # A feature amplification of a region containing a transcript
               "transcript_amplification" ]
MODERATE_IMPACT = [ "inframe_insertion",  # An inframe non synonymous variant that inserts bases into in the coding sequence    SO:0001821    Inframe insertion    MODERATE
                    "inframe_deletion",  # An inframe non synonymous variant that deletes bases from the coding sequence    SO:0001822    Inframe deletion    MODERATE
                    "missense_variant",  # A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved    SO:0001583    Missense variant    MODERATE
                    "protein_altering_variant"]
LOW_IMPACT = [
            "splice_region_variant",  #    A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron    SO:0001630    Splice region variant    LOW
            "incomplete_terminal_codon_variant",  #    A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed    SO:0001626    Incomplete terminal codon variant    LOW
            "start_retained_variant",  #    A sequence variant where at least one base in the start codon is changed, but the start remains    SO:0002019    Start retained variant    LOW
            "stop_retained_variant",  #    A sequence variant where at least one base in the terminator codon is changed, but the terminator remains    SO:0001567    Stop retained variant    LOW
            "synonymous_variant"]  #    A sequence variant where there is no resulting change to the encoded amino acid    SO:0001819    Synonymous variant    LOW
MODIFIER_IMPACT = [ "coding_sequence_variant",  #    A sequence variant that changes the coding sequence    SO:0001580    Coding sequence variant    MODIFIER
            "mature_miRNA_variant",  #    A transcript variant located with the sequence of the mature miRNA    SO:0001620    Mature miRNA variant    MODIFIER
            "5_prime_UTR_variant",  #    A UTR variant of the 5' UTR    SO:0001623    5 prime UTR variant    MODIFIER
            "3_prime_UTR_variant",  #    A UTR variant of the 3' UTR    SO:0001624    3 prime UTR variant    MODIFIER
            "non_coding_transcript_exon_variant",  #    A sequence variant that changes non-coding exon sequence in a non-coding transcript    SO:0001792    Non coding transcript exon variant    MODIFIER
            "intron_variant",  #    A transcript variant occurring within an intron    SO:0001627    Intron variant    MODIFIER
            "NMD_transcript_variant",  #    A variant in a transcript that is the target of NMD    SO:0001621    NMD transcript variant    MODIFIER
            "non_coding_transcript_variant",  #    A transcript variant of a non coding RNA gene    SO:0001619    Non coding transcript variant    MODIFIER
            "upstream_gene_variant",  #    A sequence variant located 5' of a gene    SO:0001631    Upstream gene variant    MODIFIER
            "downstream_gene_variant",  #    A sequence variant located 3' of a gene    SO:0001632    Downstream gene variant    MODIFIER
            "TFBS_ablation",  #    A feature ablation whereby the deleted region includes a transcription factor binding site    SO:0001895    TFBS ablation    MODIFIER
            "TFBS_amplification",  #    A feature amplification of a region containing a transcription factor binding site    SO:0001892    TFBS amplification    MODIFIER
            "TF_binding_site_variant",  #    A sequence variant located within a transcription factor binding site    SO:0001782    TF binding site variant    MODIFIER
            "regulatory_region_ablation",  #    A feature ablation whereby the deleted region includes a regulatory region    SO:0001894    Regulatory region ablation    MODERATE
            "regulatory_region_amplification",  #    A feature amplification of a region containing a regulatory region    SO:0001891    Regulatory region amplification    MODIFIER
            "feature_elongation",  #    A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence    SO:0001907    Feature elongation    MODIFIER
            "regulatory_region_variant",  #    A sequence variant located within a regulatory region    SO:0001566    Regulatory region variant    MODIFIER
            "feature_truncation",  #    A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence    SO:0001906    Feature truncation    MODIFIER
            "intergenic_variant",  #    A sequ
            ]
ALL_IMPACT = HIGH_IMPACT + MODERATE_IMPACT + LOW_IMPACT + MODIFIER_IMPACT
SEVERE_IMPACT = HIGH_IMPACT + MODERATE_IMPACT

COMPHET_AF_THRESH = 0.02
DNM_AF_THRESH = 0.001
DNM_GQ_THRESH = 30




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
#    Family
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


def serialize(rec):
    return rec.serialize()

def parse_max_af(V2AF):
    try:
        # example: global_ExAC_AF=0.318,global_gnomAD_AF=0.3779,global_WGS500_AF=0.3157,global_A1000G_AF=0.3173,global_UK10K_AF=0.3383
        af=[]
        for x in V2AF.split(","):
            try:
                val=float(x.split("=")[1])
                if not math.isnan(val):
                    af+=[val]
            except ValueError:
                # cannot parse
                continue
        if not af:
            return None
        return (max(af))
    except IndexError:
        print("Cannot parse AF value from %s. Setting AF to 1" % V2AF)
        return 1.0
    
class CandidateHit():

    def __init__(self, id, type):
        super().__init__()
        self.id = id
        self.type = type
        self.GT = {}
        self.GQ = {}
        self.CHROM = None
        self.POS = None
        self.REF = None
        self.ALT = None
        self.CSQ = None
        self.AF = None

    def from_tsv(self, row, ids):
        self.CHROM = row['CHROM']
        self.POS = row['POS']
        self.REF = row['REF']
        self.ALT = row['ALT']
        self.CSQ = row['Effects'].split(',')
        self.AF = parse_max_af(row['maxPopAF(global)'])
        for id in ids:
            if 'GT_' + id in row:
                self.GT[id] = Genotype(row['GT_' + id])
            else:
                self.GT[id] = None
            if 'GQ_' + id in row:
                self.GQ[id] = row['GQ_' + id]
            else:
                self.GQ[id] = None
        return(self)

    def from_record(self, rec, ids):
        self.CHROM = rec.CHROM
        self.POS = rec.POS
        self.REF = rec.REF
        self.ALT = ','.join(map(serialize, rec.ALT))
        self.AF = None  # TBD: get AF from VCF
        self.CSQ = None  # TBD: guess consequence from VEP field?
        for id in ids:
            self.GT[id] = Genotype(rec.call_for_sample[id].data["GT"])
            self.GQ[id] = rec.call_for_sample[id].data["GQ"]
        return(self)

    def to_record(self, samples, INFO):
        CALLS = []
        for id in samples:
            data = {}
            if id in self.GT:
                data["GT"] = self.GT[id].str
                data["GQ"] = str(self.GQ[id])
            else:
                data["GT"] = "./."  # id not in this pedigree
                data["GQ"] = "0"  # id not in this pedigree
            CALLS.append(vcfpy.Call(id, data, site=None))
        ALT = [vcfpy.parser.process_alt(None, self.REF, self.ALT)]
        return vcfpy.Record(self.CHROM, self.POS, ".", self.REF, ALT, None, ["PASS"], INFO, FORMAT=["GT", "GQ"], calls=CALLS)

    def getSamples(self):
        return ",".join(self.GT.keys())

    def getGTs(self):
        return ",".join(str(x) for x in self.GT.values())

    def getGQs(self):
        return ",".join(str(x) for x in self.GQ.values())

    def getMaxCSQ(self):
        if self.CSQ is None:
            return (None)
        maxidx = None
        for c in self.CSQ:
            if c in ALL_IMPACT:
                idx = ALL_IMPACT.index(c)
                if maxidx is None or idx < maxidx:
                    maxidx = idx
        if maxidx is None:
            return(None)
        return ALL_IMPACT[maxidx]

    def get_key(self):
        return(self.CHROM + ":" + str(self.POS) + "_" + self.REF + ">" + self.ALT)

    def __repr__(self):
        return self.__str__()  

    def __str__(self):
        ret = "\t".join([self.CHROM, str(self.POS), self.REF, self.ALT, self.type])
        for k in self.GT.keys():
            ret += "\t" + k + "=" + str(self.GT[k])
        return (ret)

    
class CandidateGene():

    def __init__(self, PED, GENE):
        self.PED = PED
        self.GENE = GENE
        self.CHROM = None
        self.ids = []
        self.candidates = []

    def add_candidates(self, h1, h2):
        self.candidates.append(h1)
        self.candidates.append(h2)
        self.ids.append("[" + str(h1.id) + "," + str(h2.id) + "]")
        if self.CHROM is None:
            self.CHROM = h1.CHROM
        if (self.CHROM != h1.CHROM) or (self.CHROM != h2.CHROM):
            sys.exit("Chrom mismatch! %s : %s | %s | %s " % (self.GENE, self.CHROM, h1.CHROM, h2.CHROM))

    def __repr__(self):
        return self.__str__()  

    def header():
        return ("Pedigree\tGene\tChr\tCandidateHits\tids")

    def __str__(self):
        ret = "\t".join([self.PED, self.GENE, self.CHROM, str(len(self.candidates)), ", ".join(str(x) for x in self.ids) ])
        return (ret)        
        
# print(Genotype("1/1"))
# print(Genotype(".|1"))
# print(Genotype("./."))
# print(Genotype("1/2"))
# print(Genotype("2/1"))
# print(Genotype("1/1"))
# print(Genotype("2/2"))
# print(Genotype("0/1"))
# print(Genotype("1/0"))
# print(Genotype("0/0"))
# exit(0)

def check_dnm_hit(gt, gq, gq_thresh):
    if gt is None:
        return False
    if not gt.is_hom() and not gt.is_het():
        return False
    if gq is not None and gq < gq_thresh:
        return False
    return True
def check_dnm_no_hit(gt, gq, gq_thresh):
    if gt is not None and gt not in [TYPE.HOMREF, TYPE.NOCALL]:
        return False
    if gq is not None and gq < gq_thresh:
        return False
    return True    
    

#============================================================================
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-p", "--ped", type=existing_file, required=True, nargs='+', dest="peds", metavar="peds", help="Input PED file(s)")
parser.add_argument("-s", "--snv", type=existing_file, required=True, dest="snvF", metavar="snvF", help="Input SNV VCF or pre-filtered V2 TSV file")
parser.add_argument("-c", "--cnv", type=existing_file, required=True, dest="cnvF", metavar="cnvF", help="Input CNV VCF file")
parser.add_argument("-g", "--gff", type=existing_file, required=True, dest="gffF", metavar="gffF", help="Input GFF3 annotation file (used to filter intronic variants)")
parser.add_argument("-o", "--out", type=str, required=True, dest="outdir", metavar="outdir", help="output folder")
parser.add_argument("-chr", "--chrom", type=str, required=False, dest="chrom", help="Optional comma-separated list of chroms. If set, only those chroms will be processed")

parser.add_argument("-af_chv", "--comphet_af_thresh", type=float, required=False, dest="comphet_af_thresh", default=COMPHET_AF_THRESH, help="Threshold for filtering comphet candidate hits based on maximum population allele frequency")
parser.add_argument("-af_dnm", "--dnm_af_thresh",     type=float, required=False, dest="dnm_af_thresh",     default=DNM_AF_THRESH,     help="Threshold for filtering dnm candidate hits based on maximum population allele frequency")
parser.add_argument("-gq_dnm", "--dnm_gq_thresh",     type=int, required=False,   dest="dnm_gq_thresh",     default=DNM_GQ_THRESH,     help="Threshold for filtering dnm candidate hits based on maximum genotype quality")


# parser.add_argument("-t", "--maxthreads", type=str, dest="maxthreads", default=DEF_MAX_THREADS, help="maximum threads [default: %(default)s]")
# parser.add_argument("-c", "--config", type=existing_file, required=True, dest="config", help="Configuration file [default: %(default)s]")
# 
args = parser.parse_args() 
#============================================================================

chroms = None
if args.chrom:
    chroms = args.chrom.split(",")
    print("Handling chromosomes %s" % (", ".join(chroms)))
# with vcfpy.Reader.from_path(args.snvF) as reader:
#     with vcfpy.Writer.from_path('/dev/stdout', reader.header) as writer:
#         ids = reader.header.samples.names
#         for r1 in reader:
#             print(r1)
#             writer.write_record(r1)
#             INFO = {}
#             INFO["a"]="b"
#             r2 = CandidateHit(id,"type").from_record(r1, ids).to_record(INFO)
#             id+=1
#             writer.write_record(r2)
#             break
# #    print(CandidateHit("type").from_tsv(""))
# exit(1)

# ensure dirs
if not os.path.exists(args.outdir):
        print("Creating dir " + args.outdir)
        os.makedirs(args.outdir)
log = os.path.join(args.outdir, pipename.replace(" ", "-") + ".log")
if not args.outdir.endswith("/"):
    args.outdir += "/"

# read pedigrees
peds = []
for f in args.peds:
    peds += [Pedigree(f)]
sid2ped = {}
for ped in peds:
    for fam in ped.families.values():
        for p in fam.members:
            sid2ped[p] = ped
found_peds = set()  # for storing pedigrees with at least one sample in CNV or SNV dataset.

# read all exons from GFF
gff = pr.read_gff3(args.gffF)
exons = gff[gff.Feature == 'exon']
print("Loaded %i exon annotations from %s" % (len(exons), args.gffF))

# prepare hit dics
hitsFromMum = {}
hitsFromDad = {}
denovoHits = {}
for ped in peds:
    hitsFromMum[ped.name] = collections.defaultdict(list)
    hitsFromDad[ped.name] = collections.defaultdict(list)
    denovoHits[ped.name] = collections.defaultdict(list)
chunksize = 1000000
id = 1
did = 1
collected_peds = False
# =================================================================================
# iterate SNVs and collect inherited hits
# =================================================================================
if args.snvF.endswith("vcf.gz"):
    # =================================================================================
    # parse directly from VCF (slow)
    # =================================================================================
    nline = sum(1 for i in gzip.open(args.snvF, 'rb'))
    print("Parsing ~%i SNVs from %s" % (nline, args.snvF))
    # read VCF header via vcfpy (# full file: 13,079,823 SNVs)
    reader = vcfpy.Reader.from_path(args.snvF)
    hsize = len(reader.header.lines) + 1
    parser = vcfpy.parser.RecordParser(reader.header, reader.header.samples)
    
    # collect used pedigrees
    for s in reader.header.samples.names:
        if s in sid2ped:
            found_peds.add(sid2ped[s])
     
    names = ['CHROM', 'POS', 'ID', 'REF', 'ALT',
               'QUAL', 'FILTER', 'INFO', 'FORMAT'] 
    dtype = {'CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
       'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT':str}
    for s in reader.header.samples.names:
        names = names + [s]
        dtype[s] = str
    bar = progressbar.ProgressBar(max_value=nline)     
    snvs = 0  
    for d in pd.read_csv(args.snvF,
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
        # pass only
        d = d.query("FILTER=='PASS' | FILTER=='.'")
        if chroms:
            d = d[d['CHROM'].isin(chroms)]
         
        for index, row in d.iterrows():
            snvs += 1
            bar.update(snvs)
            # parse VCF record
            rec = parser.parse_line("\t".join(map(str, row)))
             
            # get overlapping genes
            overlapping_genes = []
            # filter based on consequence
            if "ANN" in rec.INFO:
                for a in rec.INFO["ANN"]:
                    for aa in rec.ALT:
                        dat = a.split("|")
                        # print(dat[0]+"/"+aa.serialize() + "-" + dat[1])
                        if dat[0] == aa.serialize() and dat[1] in SEVERE_IMPACT:
                            overlapping_genes += [dat[3]]
            if len(overlapping_genes) > 0:
                # print("Overlapping genes: %s" % (", ".join(overlapping_genes)))
                for ped in peds:
                    for fam in ped.families.values():
                        # check all affected members and count candidates from mum or dad
                        cand_mum = 0
                        cand_dad = 0
                        cand_dnm = 0
                        for p in fam.affected_members():
                            gt_c = Genotype(rec.call_for_sample[p.id].data["GT"])
                            gt_m = Genotype(rec.call_for_sample[p.mum.id].data["GT"]) if p.mum else None
                            gt_d = Genotype(rec.call_for_sample[p.dad.id].data["GT"]) if p.dad else None
                            gq_c = rec.call_for_sample[p.id].data["GQ"]
                            gq_m = rec.call_for_sample[p.mum.id].data["GQ"] if p.mum else None
                            gq_d = rec.call_for_sample[p.dad.id].data["GQ"] if p.dad else None                        
    
                            if gt_c.is_het() and gt_m.is_het() and gt_d.type() in [TYPE.HOMREF, TYPE.NOCALL]:
                                cand_mum = cand_mum + 1
                            if gt_c.is_het() and gt_d.is_het() and gt_m.type() in [TYPE.HOMREF, TYPE.NOCALL]:
                                cand_dad = cand_dad + 1
                            # add DNM candidates if not inherited and hight quality calls
                            if check_dnm_hit(gt_c, gq_c, args.dnm_gq_thresh) and check_dnm_no_hit(gt_m, gq_m, args.dnm_gq_thresh) and check_dnm_no_hit(gt_d, gq_d, args.dnm_gq_thresh):
                                cand_dnm = cand_dnm + 1
                        if cand_mum == len(fam.affected_members()):  # all affected samples contain a HET hit inherited from their mum
                            for g in overlapping_genes:
                                hitsFromMum[ped.name][g].append(CandidateHit("m"+id, "SNV" if rec.is_snv() else "INDEL").from_record(rec, fam.members.keys()))
                                id += 1
                        if cand_dad == len(fam.affected_members()):  # all affected samples contain a HET hit inherited from their dad
                            for g in overlapping_genes:
                                hitsFromDad[ped.name][g].append(CandidateHit("d"+id, "SNV" if rec.is_snv() else "INDEL").from_record(rec, fam.members.keys()))  
                                id += 1
                        if cand_dnm:
                            isDnm = True
                            for p in fam.unaffected_members():
                                # add DNM candidates based on min GQ in unaffected samples!
                                gt = Genotype(rec.call_for_sample[p.id].data["GT"]) if p.id in samples else None
                                gq = rec.call_for_sample[p.id].data["GQ"] if p.id in samples else None
                                if not check_dnm_no_hit(gt, gq, args.dnm_gq_thresh):
                                    isDnm = False
                            if isDnm:
                                denovoHits[ped.name][g].append(CandidateHit("dnm"+did, "SNV" if rec.is_snv() else "INDEL").from_record(rec, fam.members.keys()))  
                                did += 1
                # print("Genes with SNV candidate hits from mum/dad: %s/%s" % (len(hitsFromMum[ped.name]), len(hitsFromDad[ped.name])))                 
        snvs = min(bar.max_value, snvs + (chunksize - len(d.index)))
        bar.update(snvs)
elif args.snvF.endswith("tsv.gz"):
    # =================================================================================
    # parse from V2 TSV (fast)
    # =================================================================================
    nline = sum(1 for i in gzip.open(args.snvF, 'rb'))
    print("Parsing ~%i SNVs from %s" % (nline, args.snvF))   
    bar = progressbar.ProgressBar(max_value=nline)     
    snvs = 0  
    for d in pd.read_csv(args.snvF,
                         delimiter='\t',
                         encoding='utf-8',
                         header=0,
                         comment="#",
                         float_precision='high',
                         chunksize=chunksize,
                         quoting=csv.QUOTE_NONE,
                         error_bad_lines=False): 
        # pass only
        d = d.query("FILTER=='PASS' & (maxImpact=='MODERATE' | maxImpact=='HIGH')")
        if chroms:
            d = d[d['CHROM'].isin(chroms)]

        for index, row in d.iterrows():
            # collect used pedigrees
            if not collected_peds:
                collected_peds = True
                for col in list(d.columns.values):
                    if col.startswith('GT_'):
                        if col[3:] in sid2ped:
                            found_peds.add(sid2ped[col[3:]])
            snvs += 1
            bar.update(snvs)
                        
            # get overlapping genes
            overlapping_genes = row["Genes"].split(",")
            if len(overlapping_genes) > 0:
                # print("Overlapping genes: %s" % (", ".join(overlapping_genes)))
                for ped in peds:
                    for fam in ped.families.values():
                        #print( "checking  %i genes in fam %s" % (len(overlapping_genes), fam.id))

                        # check all affected members and count candidates from mum or data
                        cand_mum = 0
                        cand_dad = 0
                        cand_dnm = 0
                        for p in fam.affected_members():
                            gt_c = Genotype(row["GT_" + p.id]) if "GT_" + p.id in row else None
                            gt_m = Genotype(row["GT_" + p.mum.id]) if p.mum and "GT_" + p.mum.id in row else None
                            gt_d = Genotype(row["GT_" + p.dad.id]) if p.dad and "GT_" + p.dad.id in row else None
                            gq_c = row["GQ_" + p.id] if "GQ_" + p.id in row else None
                            gq_m = row["GQ_" + p.mum.id] if p.mum and "GQ_" + p.mum.id in row else None
                            gq_d = row["GQ_" + p.dad.id] if p.dad and "GQ_" + p.dad.id in row else None                        

                            if gt_c is not None and gt_m is not None:
                                if gt_c.is_het() and gt_m.is_het() and (gt_d is None or gt_d.type() in [TYPE.HOMREF, TYPE.NOCALL]):
                                    cand_mum = cand_mum + 1
                            if gt_c is not None and gt_d is not None:       
                                if gt_c.is_het() and gt_d.is_het() and (gt_m is None or gt_m.type() in [TYPE.HOMREF, TYPE.NOCALL]):
                                    cand_dad = cand_dad + 1
                            # add DNM candidates if not inherited and hight quality calls
                            if check_dnm_hit(gt_c, gq_c, args.dnm_gq_thresh) and check_dnm_no_hit(gt_m, gq_m, args.dnm_gq_thresh) and check_dnm_no_hit(gt_d, gq_d, args.dnm_gq_thresh):
                                cand_dnm = cand_dnm + 1

                        if cand_mum == len(fam.affected_members()):  # all affected samples contain a HET hit inherited from their mum
                            for g in overlapping_genes:
                                # create candidate from row
                                hit = CandidateHit(id, row["TYPE"]).from_tsv(row, fam.members)
                                isCand=True
                                # filter if we know the AF and it is above a threshold
                                if hit.AF is not None and hit.AF > args.comphet_af_thresh:
                                    isCand = False
                                if hit.getMaxCSQ() not in SEVERE_IMPACT:
                                    isCand = False
                                    print("not severe dnm: '%s'" % hit.getMaxCSQ())
                                # filter non-severe consequences
                                if isCand:
                                    hitsFromMum[ped.name][g].append(hit)
                                    id += 1
                        if cand_dad == len(fam.affected_members()):  # all affected samples contain a HET hit inherited from their dad
                            for g in overlapping_genes:
                                # create candidate from row
                                hit = CandidateHit(id, row["TYPE"]).from_tsv(row, fam.members)
                                isCand=True
                                # filter if we know the AF and it is above a threshold
                                if hit.AF is not None and hit.AF > args.comphet_af_thresh:
                                    isCand = False
                                if hit.getMaxCSQ() not in SEVERE_IMPACT:
                                    isCand = False
                                    print("not severe dnm: '%s'" % hit.getMaxCSQ())
                                # filter non-severe consequences
                                if isCand:
                                    hitsFromDad[ped.name][g].append(hit)
                                    id += 1
                        if cand_dnm:
                            isDnm = True
                            for p in fam.unaffected_members():
                                gt = Genotype(row["GT_" + p.id]) if "GT_" + p.id in row else None
                                gq = row["GQ_" + p.id] if "GQ_" + p.id in row else None
                                if not check_dnm_no_hit(gt, gq, args.dnm_gq_thresh):
                                    isDnm = False
                            # create candidate from row
                            hit = CandidateHit(id, row["TYPE"]).from_tsv(row, fam.members)
                            # filter non-severe consequences
                            if hit.getMaxCSQ() not in SEVERE_IMPACT:
                                isDnm = False
                                print("not severe dnm: '%s'" % hit.getMaxCSQ())
                             # filter if we know the AF and it is above a threshold
                            if hit.AF is not None and hit.AF > args.dnm_af_thresh:
                                isDnm = False
                            if isDnm:
                                denovoHits[ped.name][g].append(hit)  
                                did += 1

                # print("Genes with SNV candidate hits from mum/dad: %s/%s" % (len(hitsFromMum[ped.name]), len(hitsFromDad[ped.name])))                 
        snvs = min(bar.max_value, snvs + (chunksize - len(d.index)))
        bar.update(snvs)
else:
    sys.exit("Unknown file type for SNVs: " + args.snvF)

# =================================================================================
# iterate CNVs and add inherited hits
# =================================================================================
nline = sum(1 for i in gzip.open(args.cnvF, 'rb'))

print("Parsing ~%i CNVs from %s" % (nline, args.cnvF))
# read VCF header via vcfpy (# full file: 13,079,823 SNVs)
reader = vcfpy.Reader.from_path(args.cnvF)
# fix issue with missing CN header in CNV vcf
reader.header.add_format_line(vcfpy.OrderedDict([('ID', 'CN'), ('Description', 'Copynumber'), ('Number', '1'), ('Type', 'Float')]))
hsize = len(reader.header.lines) + 1
parser = vcfpy.parser.RecordParser(reader.header, reader.header.samples)

names = ['CHROM', 'POS', 'ID', 'REF', 'ALT',
           'QUAL', 'FILTER', 'INFO', 'FORMAT'] 
dtype = {'CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
   'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT':str}
for s in reader.header.samples.names:
    names = names + [s]
    dtype[s] = str
bar = progressbar.ProgressBar(max_value=nline)     
cnvs = 0 
samples = reader.header.samples.names
# collect used pedigrees
for s in reader.header.samples.names:
    if s in sid2ped:
        found_peds.add(sid2ped[s])
for d in pd.read_csv(args.cnvF,
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
    # pass only
    d = d.query("FILTER=='PASS' | FILTER=='.'")
    if chroms:
        d = d[d['CHROM'].isin(chroms)]
    for index, row in d.iterrows():
        cnvs += 1
        bar.update(cnvs)
        # parse VCF record
        rec = parser.parse_line("\t".join(map(str, row)))
        type = None
        for a in rec.ALT:
            if a.type == "BND":  # skip general breakpoint annotations as too unspecific 
                 continue
            type = a.value if a.value != "None" else a.type
        if type is None:
            # ignore if type cannot be identified.
            continue

        # query overlapping exons and get list of gene ids
        overlapping_genes = []
        if ("POS" in rec.INFO) and ("END" in rec.INFO):
            pos = int(rec.INFO.get("POS"))
            end = int(rec.INFO.get("END"))
            test = pr.PyRanges(chromosomes=[rec.CHROM], starts=[min(pos, end)], ends=[max(pos, end)])
            join = test.join(exons)
            if 'gene_name' in join.columns:
                overlapping_genes = set(join.gene_name.values)
#         if "Gene_name" in rec.INFO:
#             gstr = rec.INFO.get("Gene_name")
#             if gstr is not "0" and gstr is not ".":
#                 overlapping_genes=gstr.split("/")
        
        # print("Overlapping genes: %s" % (", ".join(overlapping_genes)))
        if len(overlapping_genes) > 0:
            for ped in peds:
                for fam in ped.families.values():
                    # check all affected members and count candidates from mum or data
                    cand_mum = 0
                    cand_dad = 0
                    cand_dnm = 0
                    for p in fam.affected_members():
                        gt_c = Genotype(rec.call_for_sample[p.id].data["GT"]) if p.id in samples else None
                        gt_m = Genotype(rec.call_for_sample[p.mum.id].data["GT"]) if p.mum and p.mum.id in samples else None
                        gt_d = Genotype(rec.call_for_sample[p.dad.id].data["GT"]) if p.dad and p.dad.id in samples else None
                        # print("%s %s %s" %(gt_c, gt_m, gt_d ) ) 
                        if gt_c is not None and gt_m is not None:
                            if gt_c.is_het() and gt_m.is_het() and (gt_d is None or gt_d.type() in [TYPE.HOMREF, TYPE.NOCALL]):
                                cand_mum = cand_mum + 1
                        if gt_c is not None and gt_d is not None:
                            if gt_c.is_het() and gt_d.is_het() and (gt_m is None or gt_m.type() in [TYPE.HOMREF, TYPE.NOCALL]):
                                cand_dad = cand_dad + 1
                         # add DNM candidates based on min GQ in unaffected samples!
                        if gt_c is not None and (gt_c.is_hom() or gt_c.is_het()):
                            cand_dnm = cand_dnm + 1

                    # print("mum: %i, dad %i" %(cand_mum, cand_dad ) )    
                    if cand_mum == len(fam.affected_members()):  # all affected samples contain a HET hit inherited from their mum
                        for g in overlapping_genes:
                            hitsFromMum[ped.name][g].append(CandidateHit(id, type).from_record(rec, fam.members))
                            id += 1
                    if cand_dad == len(fam.affected_members()):  # all affected samples contain a HET hit inherited from their dad
                        for g in overlapping_genes:
                            hitsFromDad[ped.name][g].append(CandidateHit(id, type).from_record(rec, fam.members))
                            id += 1
                    if cand_dnm:
                        isDnm = True
                        for p in fam.unaffected_members():
                            gt = Genotype(rec.call_for_sample[p.id].data["GT"]) if p.id in samples else None
                            gq = rec.call_for_sample[p.id].data["GQ"] if p.id in samples else None
                            if gt is None or gq is None or gq < args.dnm_gq_thresh or gt.is_het() or gt.is_hom():
                                isDnm = False
                        if isDnm:
                            denovoHits[ped.name][g].append(CandidateHit(did, "SNV" if rec.is_snv() else "INDEL").from_record(rec, fam.members.keys()))  
                            did += 1

    cnvs = min(bar.max_value, cnvs + (chunksize - len(d.index)))
    bar.update(cnvs)

# =================================================================================
# find and write candidates
# =================================================================================
# find genes that are annotated on different chroms and add to ignore list. E.g., Metazoa_SRP is annotated on multiple chromosomes. 
filtered_genes = { "." }
for ped in peds:
    genes = set(hitsFromMum[ped.name]).union(set(hitsFromDad[ped.name]))
    for g in genes: 
        for h1 in hitsFromMum[ped.name][g]:
            for h2 in hitsFromDad[ped.name][g]: 
                if h1.CHROM != h2.CHROM:
                    filtered_genes.add(g)
print("Ignored genes: %s" % (",".join(filtered_genes)))         
 
# write candidates
candidates = {}
ped2gene2candgene = {}
ped2gene2dnmgene = {}
for ped in peds:
    # FIXME: create better header!
    header = reader.header
    header.samples.names = ped.all_ids()
    with vcfpy.Writer.from_path(args.outdir + "V2.comphet." + ped.name + ".vcf", header) as outvcf:
        for fam in ped.families.values():
            genes = set(hitsFromMum[ped.name]).union(set(hitsFromDad[ped.name]))
            genes = [x for x in genes if x not in filtered_genes]
            for g in genes:
                # print("%s: mum: %i, dad: %i" % (g, len(hitsFromMum[ped.name][g]), len(hitsFromDad[ped.name][g])), file=out)
                for h1 in hitsFromMum[ped.name][g]:
                    for h2 in hitsFromDad[ped.name][g]:     
                        # print("pair\t%s\n\t\t%s" %( h1, h2 ), file=out)
                        if h1 is not None and h2 is not None:
                            # we know that all affected samples contain h1 and h2. 
                            # Now check that no unaffected sample contains both.
                            unaffectedContainingBothHits = 0             
                            for p in fam.unaffected_members():
                                if p.id not in h1.GT:
                                    print("ERR could not find %s in h1" % p.id)
                                g1 = h1.GT[p.id].type() 
                                g2 = h2.GT[p.id].type() 
                                if g1 not in [TYPE.HOMREF, TYPE.NOCALL] and g2 not in [TYPE.HOMREF, TYPE.NOCALL]:
                                    unaffectedContainingBothHits += 1
                            # print("unaffected\t%i/%i" %(  unaffectedContainingBothHits, len(fam.unaffected_members()) ), file=out)
                            if unaffectedContainingBothHits == 0:  # TODO use threshold here?
                                # add candgene
                                candgene = None
                                if ped.name in ped2gene2candgene and g in ped2gene2candgene[ped.name]:
                                    candgene = ped2gene2candgene[ped.name][g]
                                else:
                                    candgene = CandidateGene(ped.name, g)
                                candgene.add_candidates(h1, h2)
                                if ped.name not in ped2gene2candgene:
                                    ped2gene2candgene[ped.name] = {}
                                ped2gene2candgene[ped.name][g] = candgene
                                # mark as used
                                candidates[h1.id] = h1
                                candidates[h2.id] = h2
                                # write to VCF 
                                INFO = {}
                                INFO["ped_id"] = [ped.name]
                                INFO["gene"] = [g]
                                INFO["type"] = h1.type
                                outvcf.write_record(h1.to_record(ped.all_ids(), INFO))
                                INFO["type"] = h2.type
                                outvcf.write_record(h2.to_record(ped.all_ids(), INFO))
                for h in denovoHits[ped.name][g]:
                    # write to VCF 
                    INFO = {}
                    INFO["ped_id"] = [ped.name]
                    INFO["gene"] = [g]
                    INFO["type"] = "DNM"
                    outvcf.write_record(h1.to_record(ped.all_ids(), INFO))
                    allhits = hitsFromMum[ped.name][g] + hitsFromDad[ped.name][g] + denovoHits[ped.name][g]
                    if len(allhits) >=2:
                        for x in range(0, len(allhits)-1):
                            for y in range(x+1, len(allhits)):
                                h1=allhits[x]
                                h2=allhits[y]
                                candgene = None
                                if ped.name in ped2gene2dnmgene and g in ped2gene2dnmgene[ped.name]:
                                    candgene = ped2gene2dnmgene[ped.name][g]
                                else:
                                    candgene = CandidateGene(ped.name, g)
                                candgene.add_candidates(h1, h2)
                                if ped.name not in ped2gene2dnmgene:
                                    ped2gene2dnmgene[ped.name] = {}
                                ped2gene2dnmgene[ped.name][g] = candgene

# =================================================================================
# write gene tables
# =================================================================================
out_genesF = args.outdir + "V2.comphet.genes.tsv"
print("Writing comphet result genes to %s" % out_genesF)
with open(out_genesF, 'w') as out_genes:
    print(CandidateGene.header(), file=out_genes)
    for ped in peds:
        # ---------------------
        # gene table
        # ---------------------
        if ped.name in ped2gene2candgene:
            genes = ped2gene2candgene[ped.name].keys()
            for g in genes:
                candgene = ped2gene2candgene[ped.name][g]
                print(candgene, file=out_genes)

out_genesF = args.outdir + "V2.dnm.genes.tsv"
print("Writing dnm result genes to %s" % out_genesF)
with open(out_genesF, 'w') as out_genes:
    print(CandidateGene.header(), file=out_genes)
    for ped in peds:
        # ---------------------
        # gene table
        # ---------------------
        if ped.name in ped2gene2dnmgene:
            genes = ped2gene2dnmgene[ped.name].keys()
            for g in genes:
                candgene = ped2gene2dnmgene[ped.name][g]
                print(candgene, file=out_genes)

# =================================================================================
# write hits table
# =================================================================================
outF = args.outdir + "V2.comphet.hits.tsv"
print("Writing results to %s" % outF)
with open(outF, 'w') as out:
    # ---------------------
    # individual hits 
    # ---------------------
    print("ID\tPedigree\tGene\tType\tChr\tPos\tRef\tAlt\tCSQ\tmaxAF\tGT\tGQ", file=out)
    for id in sorted(candidates.keys()):
        h1 = candidates[id]
        print("%i\t%s\t%s\t%s\t%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s" % 
              (id, ped.name, g,
               h1.type, h1.CHROM, h1.POS, h1.REF, h1.ALT,
               h1.getMaxCSQ(), str(h1.AF), h1.getGTs(), h1.getGQs()), file=out)

# =================================================================================
# find and write DNM candidates
# =================================================================================

# find genes that are annotated on different chroms and add to ignore list. E.g., Metazoa_SRP is annotated on multiple chromosomes. 
filtered_genes = { "." }
for ped in peds:
    genes = set(denovoHits[ped.name])
    for g in genes: 
        for h in denovoHits[ped.name][g]:
            for h1 in hitsFromMum[ped.name][g]: 
                for h2 in hitsFromDad[ped.name][g]: 
                    if h != h1.CHROM or h1.CHROM != h2.CHROM:
                        filtered_genes.add(g)
print("Ignored genes: %s" % (",".join(filtered_genes))) 
 # write candidates

outF = args.outdir + "V2.dnm.hits.tsv"
print("Writing DNM results to %s" % outF)
with open(outF, 'w') as out:
    ped2gene2candgene = {}
    # write header
    # for p in peds:
    #    print("# %s" % p, file=out)
    print("ID\tPedigree\tGene\tType\tChr\tPos\tRef\tAlt\tCSQ\tmaxAF\tGT\tGQ", file=out)
    for ped in peds:
        for fam in ped.families.values():
            genes = set(denovoHits[ped.name])
            genes = [x for x in genes if x not in filtered_genes]
            for g in genes:
                for h in denovoHits[ped.name][g]:
                    otherhits = len(hitsFromMum[ped.name][g]) + len(hitsFromMum[ped.name][g])
                    print("%i\t%s\t%s\t%s\t%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s" % 
                          (h.id, ped.name, g,
                           h.type, h.CHROM, h.POS, h.REF, h.ALT,
                           h.getMaxCSQ(), str(h.AF), h.getGTs(), h.getGQs()), file=out)
                    id += 1

# =================================================================================
# write all pedigrees with samples in SNV or CNV dataset
# =================================================================================

# write compiled pedigree info to fil
outF = args.outdir + "V2.peds.tsv"
print("Writing pedigrees to %s" % outF)
with open(outF, 'w') as out:
    print("Pedigree\tid\tsex\taffected\thas_parent\tis_parent", file=out)
    for ped in found_peds:
        for fam in ped.families.values():
            for p in fam.members.values():
                print(ped.name + "\t" + p.id + "\t" + ("m" if p.sex == "1" else "f" if p.sex == "2" else p.sex) + "\t" + ("1" if p.affected else "0") + "\t" + ("1" if p.has_parent() else "0") + "\t" + ("1" if p.is_parent() else "0"), file=out)
       
print("Done.")
