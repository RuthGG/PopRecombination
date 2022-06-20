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
parser.add_argument("-d", "--description", type = str,required = True, metavar = "STR", default = "", help = "run description")

args = parser.parse_args()
(lmapFile, pvalues, outDir, runName) = ( args.lmap, args.pvalues, args.out, args.description)

# %%  Test arguments
# lmapFile = "../../data/SpenceSong_hg19_recMaps_processed/CEU_recombination_map_hg19_allChr.bed"
# pvalues = '0.001,0.01,0.1,0.999,0.99,0.9'
# outDir = "../../tmp"
# runName = "test"

# %%  Set directories and log info
outDir = outDir+"/GATfiles/Annotations/"

os.makedirs(outDir, exist_ok=True)

logging.basicConfig(filename=outDir+"/log_"+runName+".txt", filemode='w',format='%(asctime)s - %(message)s', level=logging.INFO)

logging.info("Starting recMapAnnotation.py")

# %% Get data
logging.info("""Parameters

    map                     : {}
    pvalues                 : {} 
    output directory        : {}
""".format(lmapFile, pvalues, outDir)) 

logging.info( "Loading data")
# -------------------------------------

recRates = pd.read_csv(lmapFile, sep = "\t", header = 0)
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

annotations.to_csv("{}/annotations_{}.bed".format(outDir, runName), sep = "\t", header=False, index = False)

# %% End
logging.info("Finished")

