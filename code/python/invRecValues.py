#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# %% 
# Import packages
# -------------------------------

import os
from pandas import Series, DataFrame
import pandas as pd
import argparse
import math
from statistics import mean
import warnings
warnings.filterwarnings('ignore')
from datetime import datetime
import logging
import numpy as np
from scipy import stats

# %% 
# Parse arguments 
# -------------------------------

parser = argparse.ArgumentParser(description='Make annotation file with ERRs for a recombination map.')
parser.add_argument("-i", "--inversions", type=str, metavar = "FILE", required=True, help="inversions csv")
parser.add_argument("-l", "--lmap", type=str,metavar = "FILE",required=True,   help="likelihood-based map")
parser.add_argument("-o", "--out", type=str,metavar = "DIR", default = "",  help="output directory")
parser.add_argument("-d", "--description", type = str, metavar = "STR", default = "", help = "run description, recommended to start with '_'")

args = parser.parse_args()
(invFile, lmapFile, outDir, runName) = (args.inversions, args.lmap,  args.out, args.description)

# %%  Test arguments
# invFile = '../../data/InversionsAnnotation_131Inv_20211117.csv'
# lmapFile = "../../data/SpenceSong_hg19_recMaps_processed/CEU_recombination_map_hg19_allChr.bed"
# outDir = "../../tmp"
# runName = "test"

# %%  Set directories and log info
outDir = outDir+"/"+datetime.now().strftime("%Y%m%d")+"_invrecValues"+runName

if not os.path.exists(outDir):
    os.makedirs(outDir)

logging.basicConfig(filename=outDir+"/log.txt", filemode='w',format='%(asctime)s - %(message)s', level=logging.INFO)

logging.info("Starting invRecValues.py")

# %% Get data
logging.info("""Parameters

    inversions              : {}
    map                     : {}
    output directory        : {}
""".format(invFile, lmapFile, outDir)) 

logging.info( "Loading data")
# -------------------------------------

# Recombination map
recRates = pd.read_csv(lmapFile, sep = "\t", header = None)
recRates.set_axis(["Chrom", "start", "end", "recrate_base_gen"], axis = 1, inplace= True)

# Inversions
inv = pd.read_csv(invFile, sep = '\t', skiprows=1, skip_blank_lines=True)
inv = inv.iloc[:,[0,1,2,11,12,13,14]]
inv = inv.dropna()

Ori_fixed =  inv["Origin"].replace(regex = ["NAHR.*"], value = "NAHR" )
Ori_fixed = Ori_fixed.replace(regex = ["NH.*"], value = "NH" )
inv["OriginFixed"] = Ori_fixed

# inv["CenterPos"] = inv["BP1_1.1"] + ((inv["BP2_2.1"]- inv["BP1_1.1"] +1 )/2) -1 
inv["size"] = inv["BP2_2.1"] - inv["BP1_1.1"] +1

# %% A function to make mean recombination rate
def meanMap(chromosome, start, end):
    
    # Take required map
    recMapChr = recRates[recRates.Chrom == chromosome].copy()

    # Take only windows overlapping with inversion
    recMap_trimleft = recMapChr[recMapChr["end"] >= start]
    recMap_trim= recMap_trimleft[recMap_trimleft["start"] <= end]
   
    # There should be one of these MAX
    recMap_trim.loc[recMap_trim.start <= start, "start"] = start-1 # 1 to 0 based
    recMap_trim.loc[recMap_trim.end >= end, "end"] = end

    # Weigth of each window
    recMap_trim["winSize"] = recMap_trim["end"]-recMap_trim["start"] # map is 0 based I CHECKED
    recMap_trim["weight"] = recMap_trim["winSize"]/ sum(recMap_trim["winSize"])
    
    # Weighted mean
    return(sum(recMap_trim["weight"] * recMap_trim["recrate_base_gen"]))

def summaryMap(chromosome):
    # Take required map
    recMapChr = recRates[recRates.Chrom == chromosome].copy()
    # Expand each window
    listedMap = np.repeat(recMapChr.recrate_base_gen,  (recMapChr["end"]-recMapChr["start"]) )
    # Return summary
    return(listedMap.describe())

def quantileInv(invsTable, recRates):
    
    invsTable_melted = pd.melt(invsTable, id_vars = ["INV", "Origin", "Chr", "BP1_1.1", "BP1_2.1", "BP2_1.1", "BP2_2.1", "OriginFixed", "size" ])
    
    for chromosome in invsTable["Chr"].unique():
        # Take required map
        recMapChr = recRates[recRates.Chrom == chromosome]
        listedMap = np.repeat(recMapChr.recrate_base_gen,  (recMapChr["end"]-recMapChr["start"]) )
        # Make pvalue
        # invsTable_melted.loc[invsTable_melted.Chr == chromosome,'pvalue'] = [stats.percentileofscore( listedMap, value, kind='mean') for value in invsTable_melted.loc[invsTable_melted.Chr == chromosome]['value']]
        invsTable_melted.loc[invsTable_melted.Chr == chromosome,'pvalue'] = [stats.percentileofscore( listedMap, value, kind='rank') for value in invsTable_melted.loc[invsTable_melted.Chr == chromosome]['value']]
        
    return(invsTable_melted)



# %% Mean recrates
logging.info( "Calculating mean recombination rates")
# -------------------------------------

# Inside
inv['inside_WMR'] = inv.apply(lambda row: meanMap(row['Chr'],row['BP1_1.1'], row['BP2_2.1']),axis=1)

# Flanking_2side_1invsize
inv['left_invSize_WMR'] = inv.apply(lambda row: meanMap(row['Chr'],row['BP1_1.1'] - row['size']+1, row['BP1_1.1']),axis=1)
inv['right_invSize_WMR'] = inv.apply(lambda row: meanMap(row['Chr'],row['BP2_2.1'] , row['BP2_2.1'] + row['size']-1 ),axis=1)

# Flanking_2side_2invsize
inv['left_2invSize_WMR'] = inv.apply(lambda row: meanMap(row['Chr'],row['BP1_1.1'] - (row['size']*2)+1, row['BP1_1.1']),axis=1)
inv['right_2invSize_WMR']= inv.apply(lambda row: meanMap(row['Chr'],row['BP2_2.1'] , row['BP2_2.1'] + (row['size']*2)-1 ),axis=1)

# Flanking_2side_10kb
inv['left_10kb_WMR']= inv.apply(lambda row: meanMap(row['Chr'],row['BP1_1.1'] - (10000)+1, row['BP1_1.1']),axis=1)
inv['right_10kb_WMR']= inv.apply(lambda row: meanMap(row['Chr'],row['BP2_2.1'] , row['BP2_2.1'] + (10000)-1 ),axis=1)

# Flanking_2side_20kb
inv['left_20kb_WMR']= inv.apply(lambda row: meanMap(row['Chr'],row['BP1_1.1'] - (20000)+1, row['BP1_1.1']),axis=1)
inv['right_20kb_WMR']= inv.apply(lambda row: meanMap(row['Chr'],row['BP2_2.1'] , row['BP2_2.1'] + (20000)-1 ),axis=1)

# %% p-values
logging.info( "Calculating p-values for mean recombination rates")
# -------------------------------------

invPartP = quantileInv(inv.loc[(inv.Chr != "chrX") & (inv.Chr != "chrY")], recRates)

# -------------------------------------
# %% chromosome summary
# logging.info( "Calculating chromosome summary")
# -------------------------------------


# For each chromosome, distribution of recombination rates
# If we have recrate/base/generation, we can represent each value a number of = their length
# And then calculate summaries
# chromCoords = DataFrame(recRates.groupby("Chrom")["recrate_base_gen"].describe())
# chromCoords = chromCoords.reset_index(level=0)
# cc = DataFrame({"Chr" : recRates["Chrom"].unique()})
# cc['count'], cc['mean'], cc['std'], cc['min'], cc['25p'], cc['50p'], cc['75p'], cc['max'] = zip(*cc['Chr'].map(summaryMap))

# %% Write files
logging.info("Saving files")

# Write result
invPartP.to_csv( "{}/inversionRates.csv".format(outDir), index=False, sep = "\t")
# cc.to_csv( "{}/chromRates.csv".format(outDir), index=False, sep = "\t")

# %% End
logging.info("Finished")