#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# %% 
# Import packages
# -------------------------------

import os
import random
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

from joblib import Parallel, delayed
import multiprocessing


# %% 
# Parse arguments 
# -------------------------------

parser = argparse.ArgumentParser(description='Divide chromosomes effective space into windows, different methods available.')
parser.add_argument("-m", "--recmap", type=str, metavar ="FILE", required=True, help="recombination map, 0-based, right-exclusive bedfile")
parser.add_argument("-b","--band", type=str,metavar = "FILE",required=True,   help="band file to calculate centromere locations, 0-based, right-exclusive")
parser.add_argument("-f", "--fixedArmWindows", type=int,metavar = "INT", default = 0, help="number of windows to divide chromosome arms into")
parser.add_argument("-s", "--COzonesStrength", type=float,metavar = "INT", default = 0, help="between 0 and 1, persistence percentage for local minima and maxima to divide intro CO zones")
parser.add_argument("-r", "--COzonesReplicas", type=int,metavar = "INT", default = 0, help="number of replicas for map resampling to divide intro CO zones")
parser.add_argument("-o", "--out", type=str,metavar = "DIR", default = "",  help="output directory")
parser.add_argument("-d", "--description", type = str, metavar = "STR", default = "", help = "run description")

args = parser.parse_args()
(recFile, bandFile, fragCount, COstrength, COreplicas, outDir, runName) = (args.recmap, args.band, args.fixedArmWindows, args.COzonesStrength, args.COzonesReplicas, args.out, args.description)

# %% 
# Test arguments
# -------------------------------

# recFile = "../../data/SpenceSong_hg19_recMaps_processed/CEU_recombination_map_hg19_allChr.bed"
# bandFile = '../../data/cytoBand.txt'
# fragCount = 0
# COstrength = 0
# COreplicas = 80000
# outDir = "../../report/20220519_tutorialLogisticModel/data"
# runName = "avgBherer"

# Static arguments
# -------------------------------
seed = random.randrange(0,2**32-1)
np.random.seed(seed)

# Buffer for centromere
buffer=15000

# %%  
# Set directories and log info
# -------------------------------
outDir = outDir+"/divideChromosomes/"
outDir_fixed = outDir+runName+"_fixedArms_"+str(fragCount)
outDir_COzones= outDir+runName+"_COzones_"+str(COstrength)+"_"+str(COreplicas)

if not os.path.exists(outDir):
    os.makedirs(outDir)

if not os.path.exists(outDir_fixed):
    os.makedirs(outDir_fixed)

if not os.path.exists(outDir_COzones):
    os.makedirs(outDir_COzones)

logging.basicConfig(filename=outDir+"log_"+runName+"_"+str(fragCount)+"_"+str(COstrength)+"_"+str(COreplicas)+".txt", filemode='w',format='%(asctime)s - %(message)s', level=logging.INFO)

logging.info("Starting divideChromosomes.py")

# %% 
# Load and clean data
# -------------------------------
paramInfo = """Parameters

    recombination map       : {}
    band file               : {}
    numer of windows        : {} 
    crossover zone strength : {}
    crossover zone replicas : {}
    output directory        : {}
    seed:                   : {}
""".format(recFile, bandFile, fragCount,COstrength, COreplicas, outDir, seed)
logging.info(paramInfo) 

logging.info( "Loading data")
# -------------------------------

# Chromosome regions, 0-based, right-exclusive 
band = pd.read_csv(bandFile, sep = '\t', header = None )
band.set_axis(['chrom','chromStart','chromEnd','name','gieStain'], axis = 1, inplace = True)

# Recombination map, 0-based, right-exclusive
recRate = pd.read_csv(recFile, sep = "\t", header = 0)
recRate.set_axis(["Chromosome", "Start", "End", "Rate"] ,axis = 1, inplace = True)
recRate["winSize"] = recRate["End"]-recRate["Start"] 

# %% 
# Set basic chromosome limits
logging.info( "Calculating theoretical chromosome limits")
# -------------------------------

# Select centromeres End as p arm End (0-based, right-exclusive)
center = band[(band['gieStain'] == 'acen') & (band['name'].str.startswith('p') )]
center = center[['chrom', 'chromEnd']]
center.set_axis(['chrom','End'], axis = 1, inplace = True)

# Select chromosome End as q arm End (0-based, right-exclusive)
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

# %% 
# Set alternative chromosome limits according to recombination rates
logging.info( "Calculating actual chromosome limits")
# ---------------------------------------------------------------------------

# Take map limits BUT it is prepared to return other variables as well if needed
def summarizeMap(recMap):
    # Classify sections according to chromosome arm
    recMap["chromArm"] = "q"
    for c in recMap.Chromosome.unique():
        # Calculate centromere
        cen = int(chromArms[(chromArms.chrom == c) & (chromArms.chromArm == "p")].End)
        # Classify by chromosome arm
        recMap.loc[ (recMap.Chromosome == c) & (recMap.End <= cen ), "chromArm"] = "p"
        # Mark low-confidence data (windows containing centromere, which are usually significantly large)
        recMap.loc[ (recMap.Chromosome == c) & (recMap.Start < cen ) & (recMap.End >= cen ) , "chromArm"] = "cen"
        # Add buffer just in case
        cenStart=recMap.loc[(recMap.Chromosome == c) & (recMap.chromArm == "cen"), "Start"]
        cenEnd=recMap.loc[(recMap.Chromosome == c) & (recMap.chromArm == "cen"), "End"]
        if len(cenStart >0):
            recMap.loc[ (recMap.Chromosome == c) & (recMap.End < int(cenEnd)+buffer ) & (recMap.Start >= int(cenStart)-buffer), "chromArm"] = "cen"
    
    recMap["chromID"] = recMap.Chromosome + recMap.chromArm 
    # Find minimum and maximum position for each chromosome arm
    limits = DataFrame([ [min(recMap[recMap.chromID == c].Start),    max(recMap[recMap.chromID == c].End),c,    min(recMap[recMap.chromID == c].Chromosome) ] for c in recMap.chromID.unique()], columns = ["Start", "End", "chromID", "Chromosome"])
    
    return(limits) # I could also return (recMap, recMap[recMap.chromArm == "cen"])

# Join all limits in the same table
chromLimits = summarizeMap(recRate)

# Write result
# chromLimits.to_csv( "{}/mapChromosomeLimits_{}.txt".format(outDir,runName), index=False, sep = "\t")
# logging.info("  File generated: {}/mapChromosomeLimits_{}.txt """.format(outDir,runName))

# %%
# Make fixed window regions
if fragCount == 0:
    logging.info("Skipping fixed window sectioning...")
else: 
    logging.info("Dividing chromosome arms into {} sections """.format(fragCount))
# ------------------------------------------

    # Only arm dataLimits
    centromeres = chromLimits["chromID"].str.contains("cen")
    dataLimits = chromLimits[~centromeres]

    # Generate fragCount number of winRegions
    winRegions = DataFrame()

    for index, row in dataLimits.iterrows():
        winSize = math.ceil((row["End"]-row["Start"])/fragCount)
        regionStarts = range(row["Start"], row["End"], winSize) # Right coordinate might surpass arm limits
        armRegions = DataFrame( [ [x,x+winSize, row["chromID"], row["Chromosome"]] for x in regionStarts], columns=["Start", "End", "chromID", "Chromosome"])
        armRegions["winID"] = armRegions["chromID"] + [str(x) for x in range(0, fragCount, 1)]
        winRegions = winRegions.append(armRegions)

    winRegions.to_csv( "{}/windows.txt".format(outDir_fixed), index=False, sep = "\t")
    chromLimits.to_csv( "{}/workspace.txt".format(outDir_fixed), index=False, sep = "\t")
    with open("{}/log.txt".format(outDir_fixed), 'w') as f:
        print(paramInfo, file=f)

    logging.info("""    Files generated: 
        {}/windows.txt
        {}/workspace.txt
        {}/log.txt """.format(outDir_fixed, outDir_fixed,outDir_fixed))

# %%
# Make Crossover Zones
if COreplicas == 0 :
    logging.info("Skipping crossover zone sectioning...")
else: 
    logging.info("Dividing chromosomes into crossover sections with {} strength and {} replicas """.format(COstrength, COreplicas))
# ------------------------------------------
# IF ANY OF THESE FUNCTIONS IS UPDATED, REMEMBER TO UPDATE IT IN MAKEISOCHORES AS WELL!!
    def histToDensity(chroms, starts, ends, values, samples=COreplicas):
        
        #### Makes ready-to-plot density values
        # starts = histogram window start coordinates 0 -based
        # ends = histogram window end coordinates 0 - based
        # values = histogram values
        # samples = number of resamples to choose randomly from histogram
        ####################################################################3

        # Height of histogram
        h = values.to_numpy(copy =True)
        # Size-weighted height of histogram
        wh = h * (ends-starts) 

        # # Resample an uneven histogram
        # - 1 make random samples of indexes according to histogram values
        indexlist = np.random.choice(starts.index, size=samples, p=wh/wh.sum(), )

        # - 2 make random samples of coordinates within each index, all coordinates same probability
        def processIndex(i):
            return np.random.choice(range(starts[i],ends[i]), size=1)
        num_cores = multiprocessing.cpu_count()
        resamples = Parallel(n_jobs=num_cores)(delayed(processIndex)(i) for i in indexlist)
        resamples = DataFrame({"locus" : Series(np.concatenate(resamples)), "chrom":  chroms[indexlist].reset_index().Chromosome})
       
        # resamples = np.concatenate([np.random.choice(range(starts[i],ends[i]), size=1) for i in indexlist])

        # Make list of chromosomes
        mydens = DataFrame()
        # Make density for each chromosome
        for chrom in chroms.unique():

            # Make positions subset
            cidx=chroms[chroms == chrom].index
            starts_c = starts[cidx]
            ends_c = ends[cidx]
            
            # Positions for density line
            x = np.linspace(starts_c.min(), ends_c.max())

            # Find resampled kde
            rkde = sts.gaussian_kde(resamples[resamples.chrom == chrom].locus)

            # Transform kde to dataframe
            mydens = mydens.append(DataFrame({"pos" : x, "val":rkde.pdf(x), "Chromosome": [chrom]*len(x) }))

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
    def makeCOzones(InputData, colname, samples=COreplicas, strength=COstrength):
        ##### Return data density, a classified table with local minima and maxima and another table with crossover zones 
        # InputData = a table, must contain Start, End and Chromosome columns
        # + all values from the same chromosome (or else it doesn't make sense)
        # + sorted coordinates
        # colname = name of the column to analyze
        # samples = number of resamples to choose randomly from histogram
        # strength = value between 0 and 1 indicating the magnitude of persistence, which will be calculated as a portion of max histogram value
        ##########################################################3

        # Make densities
        winDensities = histToDensity(InputData["Chromosome"],InputData["Start"], InputData["End"], InputData[colname], samples)
        
        myextremes= DataFrame()
        myzones= DataFrame()

        # Local extremes
        for chrom in winDensities.Chromosome.unique():
            winDensChrom = winDensities[winDensities.Chromosome == chrom]
            idxesmin, idxesmax = localminmax(winDensChrom.valScaled, max(winDensChrom.valScaled)*strength)
            localMinima=winDensChrom.iloc[idxesmin,]
            localMinima["Type"] = "Minima"
            localMaxima=winDensChrom.iloc[idxesmax,]
            localMaxima["Type"] = "Maxima"

            extremes=localMinima.append(localMaxima)
            myextremes = myextremes.append(extremes)

            # And now crossover zones! 
            chromStart=min(chromLimits[chromLimits.Chromosome==chrom].Start)
            chromEnd=max(chromLimits[chromLimits.Chromosome==chrom].End)
            posSE = pd.concat([localMinima.pos, Series([chromStart,chromEnd ])])           
            possorted = posSE.unique()
            possorted.sort()
            zones = DataFrame({"Start": possorted[:-1], "End": possorted[1:] })

            zones["Chromosome"] = chrom
            myzones = myzones.append(zones)

        return(winDensities, myextremes, myzones)

    recRate_censored = recRate[~(recRate.chromID.str.contains("cen")) ].copy()
    recRate_censored=recRate_censored.reset_index(drop=True)
    (densities, extremes, COzones) = makeCOzones(recRate_censored, "Rate", COreplicas, COstrength)


    COzones.to_csv( "{}/windows.txt".format(outDir_COzones), index=False, sep = "\t")
    extremes.to_csv( "{}/extremes.txt".format(outDir_COzones), index=False, sep = "\t")
    densities.to_csv( "{}/densities.txt".format(outDir_COzones), index=False, sep = "\t")
    chromLimits.to_csv( "{}/workspace.txt".format(outDir_COzones), index=False, sep = "\t")
    with open("{}/log.txt".format(outDir_COzones), 'w') as f:
            print(paramInfo, file=f)
    logging.info("""    Files generated: 
        {}/windows.txt
        {}/extremes.txt
        {}/densities.txt
        {}/workspace.txt
        {}/log.txt """.format(outDir_COzones,outDir_COzones,outDir_COzones,outDir_COzones,outDir_COzones))


# %% Check with a plot

# from plotnine import ggplot,ggtitle,theme,geom_rect, aes, facet_wrap, geom_line, geom_point, scale_fill_manual, scale_color_manual

# COzones["Color"]=np.resize(["a", "b"],len(COzones))


# [
#     ggplot()+
#     geom_rect(COzones, aes(xmin = "Start", xmax = "End", fill = "Color", ymin = 0, ymax = np.inf ), alpha = 0.3)+
#     # geom_rect(recRate,aes(xmin = "Start", xmax = "End",ymin=0, ymax = "Rate" ))+
#     # geom_rect(centromeres, aes(xmin = "Start", xmax = "End", ymin = 0, ymax = np.inf), alpha = 0.7)+
#     geom_line(densities,aes(x = "pos", y = "valScaled"))+
#     geom_point(extremes, aes(x = "pos", y = "valScaled", color = "Type"))+
#     facet_wrap("Chromosome", scales = "free")+
#     # scale_fill_manual(values=["#737373", "#e1e5eb"], guide=False)+
#     scale_color_manual(values = ["#bd2b43", "#2b63bd"], guide = False)+
#     ggtitle("Crossover zones in recombination map")+
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

# %% End
logging.info("Finished")

