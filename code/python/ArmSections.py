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

parser = argparse.ArgumentParser(description='Process inversions, repeats and recombination maps and divides genome into windows.')
parser.add_argument("-i", "--inversions", type=str, metavar = "FILE", required=True, help="inversions csv")
parser.add_argument("-r", "--repeats", type=str,metavar="FILE", required=True, help="file with repeats")
parser.add_argument("-p", "--pmap", type=str, metavar ="FILE", required=True, help="pedigree-based map")
parser.add_argument("-l", "--lmap", type=str,metavar = "FILE",required=True,   help="likelihood-based map")
parser.add_argument("-b","--band", type=str,metavar = "FILE",required=True,   help="band file to calculate centromere locations")
parser.add_argument("-f", "--fragCount", type=int,metavar = "INT", default = 5, help="number of windows to divide chromosome arms into")
parser.add_argument("-o", "--out", type=str,metavar = "DIR", default = "",  help="output directory")
parser.add_argument("-d", "--description", type = str, metavar = "STR", default = "", help = "run description, recommended to start with '_'")

args = parser.parse_args()
(invFile,repFile, pmapFile, lmapFile, bandFile, fragCount, outDir, runName) = (args.inversions,args.repeats, args.pmap, args.lmap, args.band, args.fragCount, args.out, args.description)

# %%  Test arguments
# invFile = '../../data/InversionsAnnotation_131Inv_20211117.csv'
# repFile = '../../data/genomicSuperDups.txt'
# pmapFile = "../../data/Bherer_Refined_genetic_map_b37_procesed/BhererAllChroms.txt"
# lmapFile = "../../data/SpenceSong_hg19_recMaps_processed/CEU_SRR.bed"
# bandFile = '../../data/cytoBand.txt'
# fragCount = 5
# outDir = "../../tmp"
# runName = "test"

# %%  Set directories and log info
outDir = outDir+"/"+datetime.now().strftime("%Y%m%d")+"_ArmSections"+runName

if not os.path.exists(outDir):
    os.makedirs(outDir)

logging.basicConfig(filename=outDir+"/log.txt", filemode='w',format='%(asctime)s - %(message)s', level=logging.INFO)

logging.info("Starting ArmSections.py")

# %% 
# Load and clean data
logging.info("""Parameters

    inversions              : {}
    repeats                 : {}
    pedigree-based map      : {}
    likelihood-based map    : {}
    band file               : {}
    numer of windows        : {} 
    output directory        : {}
""".format(invFile, repFile, pmapFile, lmapFile, bandFile, fragCount, outDir)) 

logging.info( "Loading data")
# -------------------------------

# Chromosome regions
band = pd.read_csv(bandFile, sep = '\t', header = None )
band.set_axis(['chrom','chromStart','chromEnd','name','gieStain'], axis = 1, inplace = True)

# Recombination rates
recRateP = pd.read_csv(pmapFile, sep = "\t", header = None)
recRateP.set_axis(["Chromosome", "Start", "Rate", "cM"] ,axis = 1, inplace = True)

recRatePEnd = DataFrame()
for c in recRateP.Chromosome.unique():
    part = recRateP[recRateP.Chromosome == c]
    part["End"] = part["Start"].iloc[1:,].append( Series([0]), ignore_index = True).tolist()
    recRatePEnd = recRatePEnd.append(part)
   
recRateP = recRatePEnd[recRatePEnd.End > 0]
recRateP["winSize"] = recRateP["End"]-recRateP["Start"]+1

recRateL = pd.read_csv(lmapFile, sep = "\t", header = None)
recRateL.set_axis(["Chromosome", "Start", "End", "Rate"] ,axis = 1, inplace = True)
recRateL["winSize"] = recRateL["End"]-recRateL["Start"]+1

# Repeats (genomicSuperdups)
segDups = pd.read_csv(repFile, sep = '\t'  )
## Count repeats and intrachromosomal repeats
segDups.chromCenter = segDups["chromStart"] + ((segDups["chromEnd"]- segDups["chromStart"] +1 )/2) -1
segDups.otherCenter = segDups["otherStart"] + ((segDups["otherEnd"]- segDups["otherStart"] +1 )/2) -1

# Inversions
inv = pd.read_csv(invFile, sep = '\t', skiprows=1, skip_blank_lines=True)
inv = inv.iloc[:,[0,1,2,11,12,13,14]]
inv = inv.dropna()

Ori_fixed =  inv["Origin"].replace(regex = ["NAHR.*"], value = "NAHR" )
Ori_fixed = Ori_fixed.replace(regex = ["NH.*"], value = "NH" )
inv["OriginFixed"] = Ori_fixed

inv["CenterPos"] = inv["BP1_1.1"] + ((inv["BP2_2.1"]- inv["BP1_1.1"] +1 )/2) -1 

# %% 
# Set basic chromosome limits
logging.info( "Calculating theoretical chromosome limits")
# -------------------------------

# Select centromeres as p arm End
center = band[(band['gieStain'] == 'acen') & (band['name'].str.startswith('p') )]
center = center[['chrom', 'chromEnd']]
center.set_axis(['chrom','End'], axis = 1, inplace = True)

# Select chromosome end as q arm End
groupedBand = band['chromEnd'].groupby(band['chrom'])
end = DataFrame(groupedBand.max() )
end.reset_index(level=['chrom'], inplace = True)

# Set q arm Start
end = pd.merge(end, center, on = 'chrom')
end.set_axis(['chrom','End','Start'], axis = 1, inplace = True)

# Set p arm Start
center['Start'] = 0

# Set p and q flags
center['chromArm']='p'
end['chromArm']='q'

# Join data
chromArms = pd.concat([center, end])
chromArms["chromID"]= chromArms.chrom + chromArms.chromArm

# Write result
chromArms.to_csv( "{}/theoreticalChromosomeLimits.csv".format(outDir), index=False, sep = "\t")

# %% 
# Set alternative chromosome limits according to recombination rates
logging.info( "Calculating actual chromosome limits")
# ---------------------------------------------------------------------------

# Take map limits BUT it is prepared to return other variables as well if needed
def summarizeMap(recMap, mapName):
    # Classify sections according to chromosome arm
    recMap["chromArm"] = "q"
    for c in recMap.Chromosome.unique():
        # Calculate centromere
        cen = int(chromArms[(chromArms.chrom == c) & (chromArms.chromArm == "p")].End)
        # Classify by chromosome arm
        recMap.loc[ (recMap.Chromosome == c) & (recMap.End <= cen ), "chromArm"] = "p"
        # Mark low-confidence data (windows containing centromere, which are usually significantly large)
        recMap.loc[ (recMap.Chromosome == c) & (recMap.Start <= cen ) & (recMap.End >= cen ) , "chromArm"] = "cen"
    recMap["chromID"] = recMap.Chromosome + recMap.chromArm 
    # Find minimum and maximum position for each chromosome arm
    limits = DataFrame([ [min(recMap[recMap.chromID == c].Start),max(recMap[recMap.chromID == c].End),c] for c in recMap.chromID.unique()], columns = ["Start", "End", "chromID"])
    limits["source"] = mapName
    return(limits)# I could also return (recMap, recMap[recMap.chromArm == "cen"])

# Join all limits in the same table
chromLimits = summarizeMap(recRateP,"Bherer").append(summarizeMap(recRateL,"Spence"), ignore_index=True)

# Filter out those arms that have data for only one map
chromLimits = chromLimits.groupby("chromID").filter(lambda x: len(x) != 1)
chromLimits.reset_index(drop=True, inplace=True)

# Calculate definitive limits
maxStarts = DataFrame([chromLimits.iloc[chromLimits[chromLimits.chromID == c]["Start"].idxmax()] for c in chromLimits.chromID.unique()])
minEnds = DataFrame([chromLimits.iloc[chromLimits[chromLimits.chromID == c]["End"].idxmin()] for c in chromLimits.chromID.unique()])

dataLimits = maxStarts[["Start", "chromID", "source"]].merge(minEnds[["End", "chromID", "source"]], on = "chromID", suffixes = ["Start", "End"])

# Write result
dataLimits.to_csv( "{}/actualChromosomeLimits.csv".format(outDir), index=False, sep = "\t")

# %% 
# Calculate data for n windows # Make window limits
logging.info("Dividing chromosome arms into {} sections """.format(fragCount))
# ------------------------------------------

# Only arm dataLimits
centromeres = dataLimits["chromID"].str.contains("cen")
dataLimits = dataLimits[~centromeres]

# Generate fragCount number of winRegions
winRegions = DataFrame()

for index, row in dataLimits.iterrows():
    winSize = math.ceil((row["End"]-row["Start"]+1)/fragCount)
    regionStarts = range(row["Start"], row["End"], winSize)
    armRegions = DataFrame( [ [x,x+winSize-1, row["chromID"]] for x in regionStarts], columns=["Start", "End", "chromID"])
    armRegions["winID"] = armRegions["chromID"] + [str(x) for x in range(0, fragCount, 1)]
    winRegions = winRegions.append(armRegions)

winRegions["Chromosome"] = winRegions.chromID.replace({'p$':'','q$':''}, regex=True)


# %% 
# Calculate data for n windows # Make variables
logging.info("Calculating variables for each window")
# ------------------------------------------

winRegions.reset_index(drop=True, inplace=True)

# Physical size
winRegions["Length(bp)"] = winRegions.End - winRegions.Start +1

# Inversions, repeats, maps
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
    recP_part = recRateP[(recRateP.Chromosome == row["Chromosome"]) & (recRateP.Start >= row["Start"]) & (recRateP.Start < row["End"])]
    recL_part = recRateL[(recRateL.Chromosome == row["Chromosome"]) & (recRateL.Start >= row["Start"]) & (recRateL.Start < row["End"])]

    ## Recombination rate - Number of windows
    winRegions.loc[index,"BhererRegsNum"] = len(recP_part.index)
    winRegions.loc[index,"SpenceRegsNum"] = len(recL_part.index)

    ## Recombination rate - simple means
    winRegions.loc[index,"BhererRegsMean"] = mean(recP_part.Rate) if len(recP_part.index) > 0 else np.nan
    winRegions.loc[index,"SpenceRegsMean"] = mean(recL_part.Rate) if len(recL_part.index) > 0 else np.nan

    ## Recombination rate - ponderated means
    recP_part["weight"] = recP_part["winSize"]/ sum(recP_part["winSize"])
    winRegions.loc[index,"BhererRegsWeightMean"] = sum(recP_part["weight"] * recP_part["Rate"])

    recL_part["weight"] = recL_part["winSize"]/ sum(recL_part["winSize"])
    winRegions.loc[index,"SpenceRegsWeightMean"] = sum(recL_part["weight"] * recL_part["Rate"])

# Write result
winRegions.to_csv( "{}/windowData.csv".format(outDir), index=False, sep = "\t")

# %% End
logging.info("Finished")

