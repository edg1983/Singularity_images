#!/usr/bin/env python3
'''
Created on June, 2019

@author: niko
'''


from argparse import ArgumentParser, RawDescriptionHelpFormatter
import csv, sys, os, re, gzip, time
import pandas as pd
from pandas.core.frame import DataFrame
import numpy as np


# Necessary for including python modules from a parent directory
sys.path.append("opt/V2_core")
from core.utils import *

usage = '''python pipeline.py                             

  V2 postfilter pipeline

  Copyright 2019 Niko Popitsch. All rights reserved.
  
  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0
  
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''
start_time = time.time()
success = True

NONE_STR = "."
VALUE_SEP = ","

# overloads set __str__ method for nicer serialization
class myset(set):
    def __str__(self):
        return ",".join("{0}".format(n) for n in self)
    def __repr__(self):
        return self.__str__()  
    
# check success
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


# defines transitions and transversions
TsTvMatrix = {"AG":"Ts", "GA": "Ts", 
              "CT": "Ts", "TC": "Ts",
              "AC": "Tv", "CA": "Tv",
              "AT": "Tv", "TA": "Tv",
              "CG": "Tv", "GC": "Tv",
              "GT": "Tv", "TG": "Tv"}



# selects the maximum category value
def selectMax(df, colname, min2maxValues):
    r=list(range(0,len(min2maxValues))) # [0, 1, ..., n]
    # split by comma, replace '.' value, map to number, select max, map back to name.
    return df[colname].str.split(",", expand=True).replace(".", np.nan).replace(min2maxValues, r).max(axis=1).replace(r, min2maxValues)

#
# Main method
#
def postfilter(infile, outfile, genome, configFile, samples, overwrite=False, maxthreads=1):
    if files_exist(outfile) and not overwrite:
        print("Outputfile exists!")
        return False
        
    global config
    config = loadUserConfig(configFile)
    
    # to avoid "_csv.Error: field larger than field limit (131072)" errors
    # @see https://stackoverflow.com/questions/15063936/csv-error-field-larger-than-field-limit-131072
    maxInt = sys.maxsize
    while True:
        # decrease the maxInt value by factor 10 
        # as long as the OverflowError occurs.
        try:
            csv.field_size_limit(maxInt)
            break
        except OverflowError:
            maxInt = int(maxInt/10)
            
    # expected input dataypes (autodetect others)
    dtypes={"CHROM":"object", "POS": "int64", "REF": "object", "ALT": "object", "ID": "object", "FILTER": "object", 
            "ANN[*].GENE": "object", "ANN[*].EFFECT": "object", "ANN[*].IMPACT": "object"}
    if "roi" in config and genome in config["roi"]:
        for roi in config["roi"][genome].keys():
            dtypes[roi]="bool"
    # used for column renaming   
    sampledict_gt={}
    sampledict_gq={}
    for i,v in enumerate(samples):
        sampledict_gt[i]="GT_" + v
        sampledict_gq[i]="GQ_" + v
        dtypes[sampledict_gt[i]]="object"
        dtypes[sampledict_gq[i]]="int64"
    writeHeader=True
    mode='w'
    # get list of AF fields
    AFS={}
    AFS_clean={}
    if "af" in config and genome in config["af"]:
        for id in config["af"][genome].keys():
            fields=[]
            for a, af in config["af"][genome][id].items():  
                if af[1].count(",") == 1:  # configured AC,AN and calculated AF
                    fields.append(id+"_"+a+"_AF")
                else:
                    fields.append(id+"_"+a+"_"+af[1])
            AFS[id]=fields
    print("AFS: ")
    print(AFS)        

    converted=0
    chunksize=1000000
    for df in pd.read_csv(infile,delimiter='\t',encoding='utf-8', comment="#", float_precision='high', chunksize=chunksize, na_values=["na","NA","."], dtype=dtypes): 
        headers=df.columns.values
        
        # create output data frame
        # TODO: vectorize lambda functions to speed up
        o=DataFrame()
        o["CHROM"]=df["CHROM"]
        o["POS"]=df["POS"]
        o["REF"]=df["REF"]
        o["ALT"]=df["ALT"]
        o["FILTER"]=df["FILTER"]
        o["TYPE"] = (df["REF"].str.len()+df["ALT"].str.len()==2).map({ True: "SNV", False: "INDEL"})
        o["IsKnown"] = df["ID"].isna().replace({ True: "0", False: "1"})
        o["ID"]=df["ID"]
        try:
            o["Genes"]=np.array(list(map(myset, df["ANN[*].GENE"].astype('str').str.split(",", expand=False).values))) # get unique values
        except Exception as e: 
            print(getattr(e, 'message', repr(e)))
            pos0=df.iloc[0]["CHROM"] + ":" + str(df.iloc[0]["POS"])
            pos1=df.iloc[len(df.index)-1]["CHROM"] + ":" + str(df.iloc[len(df.index)-1]["POS"])
            print("Error converting Genes between %s and %s" % (pos0, pos1))
            o["Genes"]=None
            
        try:
            o["Effects"]=list(map(myset, df["ANN[*].EFFECT"].astype('str').str.split(",", expand=False).values))
        except Exception as e:
            print(getattr(e, 'message', repr(e)))
            pos0=df.iloc[0]["CHROM"] + ":" + str(df.iloc[0]["POS"])
            pos1=df.iloc[len(df.index)-1]["CHROM"] + ":" + str(df.iloc[len(df.index)-1]["POS"])
            print("Error converting Effects between %s and %s" % (pos0, pos1))
            o["Effects"]=None
        
        try:
            o["maxImpact"]=selectMax(df, "ANN[*].IMPACT", ["UNKNOWN","MODIFIER","LOW","MODERATE","HIGH"])
        except Exception as e:
            print(getattr(e, 'message', repr(e)))
            pos0=df.iloc[0]["CHROM"] + ":" + str(df.iloc[0]["POS"])
            pos1=df.iloc[len(df.index)-1]["CHROM"] + ":" + str(df.iloc[len(df.index)-1]["POS"])
            print("Error converting maxImpact between %s and %s" % (pos0, pos1))
            o["maxImpact"]=None
        
        # create clean AFS columns
        if "af" in config and genome in config["af"]:
            for id in config["af"][genome].keys():
                aflist=[]
                clean=[]
                for a in AFS[id]:
                    if a in df.columns and not df[a].isnull().all(): # ignore Nan columns
                        if df[a].dtype==np.float64:
                            df[a+"_clean"]=df[a]
                        else:  # split comma-separated values and get max value
                            df[a+"_clean"]=df[a].astype(str).replace('nan',np.nan).str.split(",", expand=True).replace(["."], np.nan).fillna(value=np.nan).apply(pd.to_numeric).max(axis=1)
                        clean.append(a+"_clean")
                        if not "maxPopAF_"+id in o.columns:
                            o["maxPopAF_"+id]=(a+"="+df[a].astype(str))
                        else:
                            o["maxPopAF_"+id]=o["maxPopAF_"+id]+","+(a+"="+df[a].astype(str))
                # map to categories (or NaN if no value!)
                o["maxPopAF_cat_"+id]=pd.cut( df[clean].max(axis=1), bins=[ 0, 0.0001, 0.001, 0.02, 1], right=True, labels = [ "VERY_RARE", "RARE", "COMMON", "VERY_COMMON"])

        o["TsTv"]=(df["REF"].str.upper()+df["ALT"].str.upper()).map(TsTvMatrix)
        
        # add configured output fields
        for f in config["output"]["fields"]:
            if f[0] in df:
                if f[1]=="max":
                    # split by comma, replace '.' values with nan (will be ignored) and select maximum value
                    # TODO: support other operations (min, mean, etc.)
                    o[f[0]]=df[f[0]].astype(str).str.split(",", expand=True).replace([".", "na", "NA"], [np.nan,np.nan,np.nan]).astype(float).max(axis=1)
                elif f[1]=="min":
                    # split by comma, replace '.' values with nan (will be ignored) and select minimum value
                    # TODO: support other operations (min, mean, etc.)
                    o[f[0]]=df[f[0]].astype(str).str.split(",", expand=True).replace([".", "na", "NA"], [np.nan,np.nan,np.nan]).astype(float).min(axis=1)
                else:
                    o[f[0]]=df[f[0]]

        
        # add ROIs
        if ("roi" in config and genome in config["roi"]):
            for roi in config["roi"][genome].keys():
                if roi in df:
                    o["ROI_" + roi]=df[roi].map({True: 1, False: 0})
        
        # add GT and GQ fields
        gt=pd.DataFrame( df['GEN[*].GT'].astype(str).str.split(",", expand=True) ).rename(columns=sampledict_gt)
        gq=pd.DataFrame( df['GEN[*].GQ'].astype(str).str.split(",", expand=True) ).rename(columns=sampledict_gq)
        o=pd.concat([o,gt,gq], ignore_index=False, axis=1)


        # write output
        o.to_csv(outfile, sep='\t', mode=mode, header=writeHeader, index=False, na_rep=NONE_STR)
        
        # from now on: append to file
        mode='a'
        writeHeader=False

    checksuccess("CreateFinalTSV")
    
    print("Finished.")           
    return True
                
##############################################################################################
#        Commandline
##############################################################################################
if __name__ == "__main__":   
    #============================================================================
    parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--in", type=existing_file, required=True, dest="infile", metavar="infile", help="Input TSV file (GZIPPED)")
    parser.add_argument("-o", "--out", type=str, required=True, dest="outfile", help="output TSV file")
    parser.add_argument("-s", "--samples", type=str, nargs='+', required=True, dest="samples", help="Sample ids")
    parser.add_argument("-g", "--genome", type=str, required=True, dest="genome", help="genome id")
    parser.add_argument("-c", "--config", type=existing_file, required=True, dest="config", help="Configuration file [default: %(default)s]")
    parser.add_argument("--overwrite", required=False, action='store_true', dest="overwrite", help="If set, the output files will be overwritten if existing")
    args = parser.parse_args() 
    #============================================================================
    # load config and merge with default tool config to enable system-specific program locations

    
    postfilter(args.infile, args.outfile, args.genome, args.config, args.samples, args.overwrite)
