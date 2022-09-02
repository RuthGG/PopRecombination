#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# %% 
# Import packages
# -------------------------------

import os
import weakref
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

parser = argparse.ArgumentParser(description='Fill bed file with standard measurements for analysis.')
parser.add_argument("-m", "--recmap", type=str, metavar ="FILE", required=True, help="recombination map used to calculate windows, 0-based, right-exclusive bedfile")
parser.add_argument("-w", "--windows", type=str,metavar = "FILE", default = 0, help="bedfile with window coordinates")
parser.add_argument("-o", "--out", type=str,metavar = "DIR", default = "",  help="output directory")
parser.add_argument("-d", "--description", type = str, metavar = "STR", default = "", help = "run description")

args = parser.parse_args()
(recFile, windowCoords, outDir, runName) = (args.recmap, args.windows, args.out, args.description)

# %% 
# Test arguments
# ------------------------------
# recFile = "../../data/SpenceSong_hg19_recMaps_processed/CEU_recombination_map_hg19_allChr.bed"
# windowCoords = '../../analysis/defTest_LocationPatterns/slWins/slWins_avgBherer_filtered.bed'
# outDir = "../../tmp"
# runName = "test"
# invFile = '../../data/InversionsAnnotation_131Inv_20211117.csv'
# repFile = '../../data/genomicSuperDups.txt'

# %% 
# Static arguments
# ------------------------------

invFile = 'data/InversionsAnnotation_131Inv_20211117.csv'
repFile = 'data/genomicSuperDups.txt'


# %%  
# Set directories and log info
# -------------------------------
outDir = outDir+"/fillWindowMeasurements/"
outDir_files = outDir+"/"+runName+"/"

if not os.path.exists(outDir_files):
    os.makedirs(outDir_files)

logging.basicConfig(filename=outDir+"/log_"+runName+".txt", filemode='w',format='%(asctime)s - %(message)s', level=logging.INFO)

logging.info("Starting fillWindowMeasurements.py")

# %% 
# Load and clean data
logging.info("""Parameters

    recombination map       : {}
    window coordinates      : {}
    output directory        : {}
""".format(recFile, windowCoords, outDir)) 

logging.info("""Static data

    inversion coordinates   : {}
    repeat coordinates      : {}
""".format(invFile, repFile)) 

logging.info( "Loading data")
# ------------------------------

# Recombination map, 0-based, right-exclusive
recRate = pd.read_csv(recFile, sep = "\t", header = 0)
recRate.set_axis(["Chromosome", "Start", "End", "Rate"] ,axis = 1, inplace = True)
#recRate["winSize"] = recRate["End"]-recRate["Start"] 

# Window coordinates,  0-based, right-exclusive
winRegions =  pd.read_csv(windowCoords, sep = "\t", header = 0)

# Inversions 1-based
inv = pd.read_csv(invFile, sep = '\t', skiprows=1, skip_blank_lines=True)
inv = inv.iloc[:,[0,1,2,11,12,13,14]]
inv = inv.dropna()

Ori_fixed =  inv["Origin"].replace(regex = ["NAHR.*"], value = "NAHR" )
Ori_fixed = Ori_fixed.replace(regex = ["NH.*"], value = "NH" )
inv["OriginFixed"] = Ori_fixed

# Transform inversions to 0-based
inv["Start"] = inv["BP1_1.1"]-1
inv["End"] = inv["BP2_2.1"]
inv["CenterPos"] = inv["Start"] + ((inv["End"]- inv["Start"])/2) 

# Repeats (genomicSuperdups) - Assuming 1 -based
segDups = pd.read_csv(repFile, sep = '\t'  )
## Count repeats and intrachromosomal repeats
segDups["chromCenter"] = (segDups["chromStart"]-1) + ((segDups["chromEnd"]- (segDups["chromStart"]-1) )/2) 
segDups["otherCenter"] = (segDups["otherStart"]-1) + ((segDups["otherEnd"]- (segDups["otherStart"]-1) )/2) 


# %% 
# Calculate data for each window
logging.info("Calculating variables for each window")
# ------------------------------------------

# Just in case
winRegions.reset_index(drop=True, inplace=True) 

# Physical size  
#winRegions["Length(bp)"] = winRegions.End - winRegions.Start 

# Inversions, repeats, maps - EVERYTHING IS 0-BASED
for index, row in winRegions.iterrows():
    
    # Inversions (all, NH, NAHR)
    inv_part = inv[(inv.Chr == row["Chromosome"]) & (inv.CenterPos >= row["Start"]) & (inv.CenterPos < row["End"])]
    
    winRegions.loc[index,"invCenters"] = len(inv_part.index) 
    winRegions.loc[index,"NHCenters"] = len(inv_part[inv_part.OriginFixed == "NH" ].index) 
    winRegions.loc[index,"NAHRCenters"] = len(inv_part[inv_part.OriginFixed == "NAHR" ].index) 

    # Repeats (all, and intrachromosomal)
    rep_part_A = segDups[(segDups.chrom == row["Chromosome"]) & (segDups.chromCenter >= row["Start"]) & (segDups.chromCenter < row["End"])]
    rep_part_B = segDups[(segDups.otherChrom == row["Chromosome"]) & (segDups.otherCenter >= row["Start"]) & (segDups.otherCenter < row["End"])]

    winRegions.loc[index,"allRepCounts"] = len(rep_part_A.index) + len(rep_part_B.index)
    ## Intrachromosomal repeats, each copy counted different because they can be in different regions 
    winRegions.loc[index,"intraRepCounts"] = len(rep_part_A[(rep_part_A.chrom == rep_part_A.otherChrom)].index) + len(rep_part_B[(rep_part_B.chrom == rep_part_B.otherChrom)].index)

    # Recombination rates
    recRate_part = recRate[(recRate.Chromosome == row["Chromosome"]) & (recRate.Start < row["End"]) & (recRate.End > row["Start"])]
    if len(recRate_part.index) > 0:
        recRate_part.Start[recRate_part.Start == recRate_part.Start.min()] = row.Start
        recRate_part.End[recRate_part.End == recRate_part.End.max()] = row.End
        recRate_part["winSize"] = recRate_part.End - recRate_part.Start 

        winRegions.loc[index,"WAvgRate"] = np.average (recRate_part["Rate"], weights = recRate_part["winSize"] )
        winRegions.loc[index,"maxRate"] = max(recRate_part["Rate"])
        winRegions.loc[index,"minRate"] = min(recRate_part["Rate"])
    else:
        winRegions.loc[index,"WAvgRate"] = math.nan
        winRegions.loc[index,"maxRate"] = math.nan
        winRegions.loc[index,"minRate"] = math.nan
# %% 
# Write result
winRegions.to_csv( "{}/windowData.txt".format(outDir_files), index=False, sep = "\t")
logging.info("""    Files generated: {}/windowData.txt""".format(outDir_files))