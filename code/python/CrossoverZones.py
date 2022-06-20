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

# For KDE estimation with resampling
import matplotlib.pyplot as plt
import scipy.stats as sts

# For local minima and maxima
from persistence1d import RunPersistence

# %% 
# Parse arguments 
# -------------------------------

parser = argparse.ArgumentParser(description='Process inversions, repeats and recombination maps and divides genome into crossover zones.')
parser.add_argument("-i", "--inversions", type=str, metavar = "FILE", required=True, help="inversions csv")
parser.add_argument("-r", "--repeats", type=str,metavar="FILE", required=True, help="file with repeats")
parser.add_argument("-p", "--pmap", type=str, metavar ="FILE", required=True, help="pedigree-based map")
parser.add_argument("-l", "--lmap", type=str,metavar = "FILE",required=True,   help="likelihood-based map")
parser.add_argument("-b","--band", type=str,metavar = "FILE",required=True,   help="band file to calculate centromere locations")
parser.add_argument("-s", "--strength", type=float,metavar = "FLOAT", default = 0.1, help="between 0 and 1, persistence percentage for local minima and maxima")
parser.add_argument("-o", "--out", type=str,metavar = "DIR", default = "",  help="output directory")
parser.add_argument("-d", "--description", type = str, metavar = "STR", default = "", help = "run description, recommended to start with '_'")

args = parser.parse_args()
(invFile,repFile, pmapFile, lmapFile, bandFile, strength, outDir, runName) = (args.inversions,args.repeats, args.pmap, args.lmap, args.band, args.strength, args.out, args.description)

# %%  Test arguments
# invFile = '../../data/InversionsAnnotation_131Inv_20211117.csv'
# repFile = '../../data/genomicSuperDups.txt'
# pmapFile = "../../data/Bherer_Refined_genetic_map_b37_procesed/BhererAllChroms.txt"
# lmapFile = "../../data/SpenceSong_hg19_recMaps_processed/CEU_recombination_map_hg19_allChr.bed"
# bandFile = '../../data/cytoBand.txt'
# strength = 0.1
# outDir = "../../tmp"
# runName = "test"


# %%  Set directories and log info
outDir = outDir+"/"+datetime.now().strftime("%Y%m%d")+"_CrossoverZones"+runName

if not os.path.exists(outDir):
    os.makedirs(outDir)

logging.basicConfig(filename=outDir+"/log.txt", filemode='w',format='%(asctime)s - %(message)s', level=logging.INFO)

logging.info("Starting CrossoverZones.py")

# %% 
# Load and clean data
logging.info("""Parameters

    inversions              : {}
    repeats                 : {}
    pedigree-based map      : {}
    likelihood-based map    : {}
    band file               : {}
    persistence threshold   : {} 
    output directory        : {}
""".format(invFile, repFile, pmapFile, lmapFile, bandFile, strength, outDir)) 

logging.info( "Loading data")
# -------------------------------

# Chromosome regions - 0-BASED
band = pd.read_csv(bandFile, sep = '\t', header = None )
band.set_axis(['chrom','chromStart','chromEnd','name','gieStain'], axis = 1, inplace = True)

# Recombination rates - BOTH MAPS ARE 0-BASED - USE AS-IS
recRateP = pd.read_csv(pmapFile, sep = "\t", header = None)
recRateP.set_axis(["Chromosome", "Start", "Rate", "cM"] ,axis = 1, inplace = True)

recRatePEnd = DataFrame()
for c in recRateP.Chromosome.unique():
    part = recRateP[recRateP.Chromosome == c]
    part["End"] = part["Start"].iloc[1:,].append( Series([0]), ignore_index = True).tolist()
    recRatePEnd = recRatePEnd.append(part)
   
recRateP = recRatePEnd[recRatePEnd.End > 0]
recRateP["winSize"] = recRateP["End"]-recRateP["Start"]

recRateL = pd.read_csv(lmapFile, sep = "\t", header = None)
recRateL.set_axis(["Chromosome", "Start", "End", "Rate"] ,axis = 1, inplace = True)
recRateL["winSize"] = recRateL["End"]-recRateL["Start"]

# Repeats (genomicSuperdups)
segDups = pd.read_csv(repFile, sep = '\t'  )
## Count repeats and intrachromosomal repeats
segDups.chromCenter = segDups["chromStart"] + ((segDups["chromEnd"]- segDups["chromStart"] +1 )/2) -1
segDups.otherCenter = segDups["otherStart"] + ((segDups["otherEnd"]- segDups["otherStart"] +1 )/2) -1

# Inversions - 1-BASED COORDINATES
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
# # Make sections
# # size of windows
# fragCount = 100
# logging.info("Divide chromosomes into {} sections """.format(fragCount))
# # ------------------------------------------

# # I want the whole chromosome this time
# dataLimits["Chromosome"] = dataLimits.chromID.replace(["p", "q", "cen"], "", regex = True)
# chromLimits = dataLimits.groupby("Chromosome").agg(
#     Start = pd.NamedAgg(column = "Start", aggfunc = min),
#     End = pd.NamedAgg(column = "End", aggfunc = max)
# )

# # I store centromeres for later anyway
# centromeres = dataLimits[dataLimits["chromID"].str.contains("cen")]

# # Generate fragCount number of winRegions
# winRegions = DataFrame()

# for index, row in chromLimits.iterrows():
#     winSize = math.ceil((row["End"]-row["Start"]+1)/fragCount)
#     regionStarts = range(row["Start"], row["End"], winSize)
#     armRegions = DataFrame( [ [x,x+winSize-1, index] for x in regionStarts], columns=["Start", "End", "chromID"])
#     armRegions["winID"] = armRegions["chromID"]+"_" + [str(x) for x in range(0, fragCount, 1)]
#     winRegions = winRegions.append(armRegions)

# winRegions["Chromosome"] = winRegions.chromID #I keep this just in case...

# %% 
# Make histogram
# logging.info("Calculating histogram values for recombination variables")
# # ------------------------------------------

# winRegions.reset_index(drop=True, inplace=True)

# # Physical size
# winRegions["Length(bp)"] = winRegions.End - winRegions.Start +1

# # Inversions, repeats, maps
# for index, row in winRegions.iterrows():
    
#     # Recombination rates
#     recP_part = recRateP[(recRateP.Chromosome == row["Chromosome"]) & (recRateP.Start >= row["Start"]) & (recRateP.Start < row["End"])]
#     recL_part = recRateL[(recRateL.Chromosome == row["Chromosome"]) & (recRateL.Start >= row["Start"]) & (recRateL.Start < row["End"])]

#     ## Recombination rate - ponderated means
#     recP_part["weight"] = recP_part["winSize"]/ sum(recP_part["winSize"])
#     winRegions.loc[index,"BhererRegsWeightMean"] = sum(recP_part["weight"] * recP_part["Rate"])

#     recL_part["weight"] = recL_part["winSize"]/ sum(recL_part["winSize"])
#     winRegions.loc[index,"SpenceRegsWeightMean"] = sum(recL_part["weight"] * recL_part["Rate"])

# %% 
# Extreme points from KDE estimation with resampling
# Number of samples
n = 50
logging.info("Resample map {} times and calculate local relevant points with {} strength""".format(n, strength))
# ------------------------------------------

def histToDensity(starts, ends, values, samples=n):
    #### Makes ready-to-plot density values
    # starts = histogram window start coordinates 0 -based
    # ends = histogram window end coordinates 0 - based
    # values = histogram values
    # samples = number of resamples to choose randomly from histogram
    ####################################################################3

    # Height of histogram
    h = values.to_numpy(copy =True)
    # Size-weighted height of histogram
    wh = h * (starts-ends) 

    # Positions for density line
    e = starts.append(Series(ends.iloc[-1]), ignore_index = True) 
    e = e.to_numpy(copy =True)
    x = np.linspace(e.min(), e.max())

    # # Resample an uneven histogram
    # - 1 make random samples of indexes according to histogram values
    indexlist = np.random.choice(starts.index, size=1000, p=wh/wh.sum())
    # - 2 make random samples of coordinates within each index, all coordinates same probability
    resamples = np.concatenate([np.random.choice(range(starts[i],ends[i]), size=1) for i in indexlist])

    # # resample the histogram 
    # resamples = np.random.choice((e[:-1] + e[1:])/2, size=samples, p=h/h.sum())

    # Find resampled kde
    rkde = sts.gaussian_kde(resamples)

    # Transform kde to dataframe
    mydens = DataFrame({"pos" : x, "val":rkde.pdf(x)})

    # Scale kde values to the original values
    mydens["valScaled"] = mydens.val * ( mean(values)/mean(mydens.val) )

    # Results
    return(mydens)
def localminmax(InputData, minpersistence): 
    ####### Make local minima and maxima
    # Indexes must be ordered - resorted for it to work!!      
    # This simple call is all you need to compute the extrema of the given data (InputData) and their persistence.
    ############################################################

    ExtremaAndPersistence = RunPersistence(InputData)

    #~ Keep only those extrema with a persistence larger than theshold.
    Filtered = [t for t in ExtremaAndPersistence if t[1] > minpersistence]

    #~ Sort the list of extrema by persistence.
    Sorted = sorted(Filtered, key=lambda ExtremumAndPersistence: ExtremumAndPersistence[1])

    #~ Print to console
    minima = []
    maxima = []
    for i, E in enumerate(Sorted):
        strType = "Minimum" if (i % 2 == 0) else "Maximum"
        if (strType == "Minimum"):
            minima.append(E[0])
        else:
            maxima.append(E[0])
    # print("%s at index %d with persistence %g and data value %g" % (strType, E[0], E[1], InputData[E[0]]))

    return(minima, maxima)
def makeCOzones(InputData, colname, samples, strength):
    ##### Return data density, a classified table with local minima and maxima and another table with crossover zones 
    # InputData = a table, must contain Start, End and Chromosome columns
    # + all values from the same chromosome (or else it doesn't make sense)
    # + sorted coordinates
    # colname = name of the column to analyze
    # samples = number of resamples to choose randomly from histogram
    # strength = value between 0 and 1 indicating the magnitude of persistence, which will be calculated as a portion of max histogram value
    ##########################################################3

    # Make densities
    winDensities = histToDensity(InputData["Start"], InputData["End"], InputData[colname], samples)
    winDensities["Chromosome"] = InputData.iloc[0]["Chromosome"]

    # Local extremes
    idxesmin, idxesmax = localminmax(winDensities.valScaled, max(winDensities.valScaled)*strength)
    localMinima=winDensities.iloc[idxesmin,]
    localMinima["Type"] = "Minima"
    localMaxima=winDensities.iloc[idxesmax,]
    localMaxima["Type"] = "Maxima"

    extremes=localMinima.append(localMaxima)
    extremes["Chromosome"] = InputData.iloc[0]["Chromosome"]

    # And now crossover zones! 
    possorted = np.array(localMinima.pos.sort_values())
    myzones = DataFrame({"start": possorted[:-1], "end": possorted[1:] })
    myzones["color"] =  np.resize(["a", "b"], myzones.shape[0])
    myzones["Chromosome"] = InputData.iloc[0]["Chromosome"]

    return(winDensities, extremes, myzones)


extremes= DataFrame()
COzones= DataFrame()
densities = DataFrame()

recRateL.name = "Spence"
recRateP.name = "Bherer"

for stat in [recRateL,recRateP]:
    for chrom in stat.Chromosome.unique():

        subset = stat[stat.Chromosome == chrom]
        d_stat, ex_stat, Cz_stat = makeCOzones(stat, "Rate", n, strength)
        d_stat["Stat"] = stat.name
        ex_stat["Stat"] = stat.name
        Cz_stat["Stat"] = stat.name

        extremes = extremes.append(ex_stat)
        COzones = COzones.append(Cz_stat)
        densities = densities.append(d_stat)

# %% Check with a plot

from plotnine import ggplot,ggtitle,theme,geom_rect, aes, facet_wrap, geom_line, geom_point, scale_fill_manual, scale_color_manual

# [
#     ggplot()+
#     geom_rect(COzones_BRWM, aes(xmin = "start", xmax = "end", fill = "color", ymin = 0, ymax = np.inf ), alpha = 0.3)+
#     geom_rect(winRegions,aes(xmin = "Start", xmax = "End",ymin=0, ymax = "BhererRegsWeightMean" ))+
#     geom_rect(centromeres, aes(xmin = "Start", xmax = "End", ymin = 0, ymax = np.inf), alpha = 0.7)+
#     geom_line(densities_BRWM,aes(x = "pos", y = "valScaled"))+
#     geom_point(extremes_BRWM, aes(x = "pos", y = "valScaled", color = "Type"))+
#     facet_wrap("Chromosome", scales = "free")+
#     scale_fill_manual(values=["#737373", "#e1e5eb"], guide=False)+
#     scale_color_manual(values = ["#bd2b43", "#2b63bd"], guide = False)+
#     ggtitle("Crossover zones in Bherer recombination map (pedigrees)")+
#     theme(figure_size=(16, 8))  # here you define the plot size+
    
# ]

# [
#     ggplot()+
#     geom_rect(COzones_SRWM, aes(xmin = "start", xmax = "end", fill = "color", ymin = 0, ymax = np.inf ), alpha = 0.3)+
#     geom_rect(winRegions,aes(xmin = "Start", xmax = "End",ymin=0, ymax = "SpenceRegsWeightMean" ))+
#     geom_rect(centromeres, aes(xmin = "Start", xmax = "End", ymin = 0, ymax = np.inf), alpha = 0.7)+
#     geom_line(densities_SRWM,aes(x = "pos", y = "valScaled"))+
#     geom_point(extremes_SRWM, aes(x = "pos", y = "valScaled", color = "Type"))+
#     facet_wrap("Chromosome", scales = "free")+
#     scale_fill_manual(values=["#737373", "#e1e5eb"], guide=False)+
#     scale_color_manual(values = ["#bd2b43", "#2b63bd"], guide = False)+
#     ggtitle("Crossover zones in Spence recombination map (likelihood)")+
#     theme(figure_size=(16, 8))  # here you define the plot size+
    
# ]


# %% 
# Write intermediate results about crossover zones
logging.info("Write histogram, densities, extremes")
# ------------------------------------------

winRegions.to_csv( "{}/histogram.csv".format(outDir), index=False, sep = "\t")
extremes.to_csv( "{}/extremes.csv".format(outDir), index=False, sep = "\t")
densities.to_csv( "{}/densities.csv".format(outDir), index=False, sep = "\t")

# %% 
# Calculate data for n windows # Make variables
logging.info("Calculating variables for each window")
# ------------------------------------------

# changing names for simplicity, I don't need the histogram anymore
winRegions = COzones 

winRegions.reset_index(drop=True, inplace=True)
winRegions.columns=["Start", "End", "Color", "Chromosome", "Stat"]

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


# %% 
# Write intermediate results about crossover zones
logging.info("Write crossover zone values")
# ------------------------------------------
winRegions.to_csv( "{}/COzone_windowData.csv".format(outDir), index=False, sep = "\t")


# %% End
logging.info("Finished")
