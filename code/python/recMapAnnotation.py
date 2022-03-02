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

# %% 
# Parse arguments 
# -------------------------------

parser = argparse.ArgumentParser(description='Make annotation file with ERRs for a recombination map.')
parser.add_argument("-l", "--lmap", type=str,metavar = "FILE",required=True,   help="likelihood-based map")
parser.add_argument("-p","--pvalues", type=str,metavar = "STR",required=True,   help="comma-separated list of p-values to determine ERRs")
parser.add_argument("-o", "--out", type=str,metavar = "DIR", default = "",  help="output directory")
parser.add_argument("-d", "--description", type = str, metavar = "STR", default = "", help = "run description, recommended to start with '_'")

args = parser.parse_args()
(lmapFile, pvalues, outDir, runName) = ( args.lmap, args.pvalues, args.out, args.description)

# %%  Test arguments
# lmapFile = "../../data/SpenceSong_hg19_recMaps_processed/CEU_SRR.bed"
# pvalues = '0.01,0.001,0.0001'
# outDir = "../../tmp"
# runName = "test"

# %%  Set directories and log info
outDir = outDir+"/"+datetime.now().strftime("%Y%m%d")+"_recMapAnnotation"+runName

if not os.path.exists(outDir):
    os.makedirs(outDir)

logging.basicConfig(filename=outDir+"/log.txt", filemode='w',format='%(asctime)s - %(message)s', level=logging.INFO)

logging.info("Starting recMapAnnotation.py")

# %% Get data
logging.info("""Parameters

    map                     : {}
    pvalues                 : {} 
    output directory        : {}
""".format(lmapFile, pvalues, outDir)) 

logging.info( "Loading data")
# -------------------------------------

recRates = pd.read_csv(lmapFile, sep = "\t", header = None)
recRates.set_axis(["Chrom", "start", "end", "recrate_base_gen"], axis = 1, inplace= True)

pvals = [float(x) for x in pvalues.split(",")]

# %% Make annotations file
logging.info( "Calculating ERRs for each p-value")
# -------------------------------

annotations = DataFrame()

for sig in pvals:

    threshold = np.quantile(recRates["recrate_base_gen"], sig)

    if sig > 0.5:
        #  Hotspots
        namesig = "Hotspot_"+str( sig )
        recRates[namesig] = recRates["recrate_base_gen"] > threshold
    else :
        # Coldspots
        namesig = "Coldspot_"+str(sig)
        recRates[namesig] = recRates["recrate_base_gen"] < threshold

    subannot = recRates[recRates[namesig]][["Chrom", "start", "end"]]
    subannot["trackName"] = namesig
    annotations = annotations.append(subannot, ignore_index=True)

# %% Save annotations file
logging.info( "Saving annotations file")
# -------------------------------

annotations.to_csv("{}/annotations.bed", sep = "\t", header=False, index = False)

