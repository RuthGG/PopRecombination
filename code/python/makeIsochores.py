#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# %% 
# Import packages
# -------------------------------

import os
import random
from tkinter.tix import MAX
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

parser = argparse.ArgumentParser(description='Calculate isochores from densities for a file with coordinates.')
parser.add_argument("-f", "--feature", type=str, metavar ="FILE", required=True, help="feature, 0-based, right-exclusive bedfile")
parser.add_argument("-t", "--thresholdRange", type=float,metavar = "INT", default = 0, help="between 0 and 1, persistence percentage for local minima and maxima to divide intro CO zones")
parser.add_argument("-o", "--out", type=str,metavar = "DIR", default = "",  help="output directory")
parser.add_argument("-d", "--description", type = str, metavar = "STR", default = "", help = "run description")

args = parser.parse_args()
(featFile, threshold, outDir, runName) = (args.feature, args.thresholdRange, args.out, args.description)

# %% 
# Test arguments
# -------------------------------

# featFile =  "../../data/genomicSuperDups_coords.bed"
# threshold = 0.2
# outDir = "../../tmp"
# runName = "test"

# Static arguments
# ------------------------------

bandFile = 'data/cytoBand.txt'
COstrength=0

# %%  
# Set directories and log info
# -------------------------------
outDir = outDir+"/GATfiles/Isochores/"
outDir_COzones= outDir+runName+"_COzones_"+str(threshold)

os.makedirs(outDir_COzones, exist_ok=True)

logging.basicConfig(filename=outDir+"log_"+runName+"_"+str(threshold)+".txt", filemode='w',format='%(asctime)s - %(message)s', level=logging.INFO)

logging.info("Starting makeIsochores.py")

# %% 
# Load and clean data
# -------------------------------
paramInfo = """Parameters

    feature file            : {}
    threshold range         : {}
    output directory        : {}
""".format(featFile, threshold, outDir)
logging.info(paramInfo) 

logging.info( "Loading data")
# -------------------------------

# Feature bamfile, 0-based, right-exclusive
feature = pd.read_csv(featFile, sep = "\t", header = 0)
feature.set_axis(["Chromosome", "Start", "End"] ,axis = 1, inplace = True)

# Chromosome regions, 0-based, right-exclusive 
band = pd.read_csv(bandFile, sep = '\t', header = None )
band.set_axis(['chrom','chromStart','chromEnd','name','gieStain'], axis = 1, inplace = True)


# %%
# Make Feature Zones
logging.info("Dividing chromosomes into sections with {} strength """.format(COstrength))
# ------------------------------------------

# IF ANY OF THESE FUNCTIONS IS UPDATED, REMEMBER TO UPDATE IT IN divideChromosomes AS WELL!!

def histToDensity(chroms, starts, ends):
    #### Makes ready-to-plot density values
    # starts = histogram window start coordinates 0 -based
    # ends = histogram window end coordinates 0 - based
    ####################################################################3
    
    # Find centerpoints
    centers = starts + (ends-starts)/2
    
    # Skip resampling part
    resamples = DataFrame({"locus" : centers, "chrom":  chroms})
    
    # Make list of chromosomes
    mydens = DataFrame()
    # Make density for each chromosome
    for chrom in band.chrom.unique():

        # Make positions subset
        band_c = band[band.chrom == chrom]
        
        # Positions for density line
        x = np.linspace(band_c.chromStart.min(), band_c.chromEnd.max())

        # Find resampled kde
        rkde = sts.gaussian_kde(resamples[resamples.chrom == chrom].locus)

        # Transform kde to dataframe
        mydens = mydens.append(DataFrame({"pos" : x, "val":rkde.pdf(x), "Chromosome": [chrom]*len(x) }))

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
def makeCOzones(InputData, strength=COstrength):
    ##### Return data density, a classified table with local minima and maxima and another table with crossover zones 
    # InputData = a table, must contain Start, End and Chromosome columns
    # + all values from the same chromosome (or else it doesn't make sense)
    # + sorted coordinates
    # colname = name of the column to analyze
    # samples = number of resamples to choose randomly from histogram
    # strength = value between 0 and 1 indicating the magnitude of persistence, which will be calculated as a portion of max histogram value
    ##########################################################3

    # Make densities
    winDensities = histToDensity(InputData["Chromosome"],InputData["Start"], InputData["End"])
    
    myextremes= DataFrame()
    myzones= DataFrame()

    # Local extremes
    for chrom in winDensities.Chromosome.unique():
        winDensChrom = winDensities[winDensities.Chromosome == chrom]
        idxesmin, idxesmax = localminmax(winDensChrom.val, max(winDensChrom.val)*strength)
        localMinima=winDensChrom.iloc[idxesmin,]
        localMinima["Type"] = "Minima"
        localMaxima=winDensChrom.iloc[idxesmax,]
        localMaxima["Type"] = "Maxima"

        extremes=localMinima.append(localMaxima)
        myextremes = myextremes.append(extremes)

        # And now crossover zones! 
        possorted = np.array(localMinima.pos.sort_values())
        zones = DataFrame({"Start": possorted[:-1], "End": possorted[1:] })
        zones["Chromosome"] = chrom
        myzones = myzones.append(zones)

    return(winDensities, myextremes, myzones)

(densities, extremes, COzones) = makeCOzones(feature, COstrength)

mn=min(extremes[extremes["Type"] == "Minima"].val)
mx=max(extremes[extremes["Type"] == "Maxima"].val)
turningPoint=threshold
threshold=mn+(mx-mn)*turningPoint

myzones=DataFrame()
# And now crossover zones! 
for chrom in densities.Chromosome.unique():
    densities_c = densities[densities.Chromosome == chrom]
    densities_c["thresStatus"]= ["High" if value>=threshold else "Low" for value in densities_c["val"] ]
    vec = pd.concat([Series(["Start"]),densities_c["thresStatus"][:-2], Series(["End"])], ignore_index=True)
    densities_c=densities_c.reset_index(drop=True)
    densities_c["prevStatus"]=vec

    densities_c["turning"]=densities_c["thresStatus"] != densities_c["prevStatus"]
   
    possorted=densities_c[(densities_c["turning"]== True)]["pos"].astype(int)
    types =densities_c[(densities_c["turning"]== True)]["thresStatus"]
    zones = DataFrame({"Start": possorted[:-1].reset_index(drop=True), "End": possorted[1:].reset_index(drop=True), "Color":types[:-1].reset_index(drop=True) },)
    zones["Chromosome"] = chrom
    myzones = myzones.append(zones)


myzones[["Chromosome", "Start", "End", "Color"]].to_csv( "{}/windows.txt".format(outDir_COzones), header=False, index=False, sep = "\t")
extremes.to_csv( "{}/extremes.txt".format(outDir_COzones), index=False, sep = "\t")
densities.to_csv( "{}/densities.txt".format(outDir_COzones), index=False, sep = "\t")

with open("{}/log.txt".format(outDir_COzones), 'w') as f:
        print(paramInfo, file=f)
logging.info("""    Files generated: 
    {}/windows.txt
    {}/extremes.txt
    {}/densities.txt""".format(outDir_COzones,outDir_COzones,outDir_COzones,outDir_COzones))


# %% Check with a plot

from plotnine import ggplot,ggtitle,geom_hline,theme,geom_rect, aes, facet_wrap, geom_line, geom_point, scale_fill_manual, scale_color_manual
COzones=myzones

# COzones["Color"]=np.resize(["a", "b"],len(COzones))


p=ggplot(COzones)+\
    geom_rect(COzones, aes(xmin = "Start", xmax = "End", fill = "Color", ymin = 0, ymax = np.inf ), alpha = 0.3)+    geom_line(densities,aes(x = "pos", y = "val"))+\
    geom_point(extremes, aes(x = "pos", y = "val", color = "Type"))+\
    facet_wrap("Chromosome", scales = "free")+\
    geom_hline(yintercept=threshold)+\
    scale_color_manual(values = ["#bd2b43", "#2b63bd"], guide = False)+\
    ggtitle("Isochores for "+runName+" with threshold = {} ({})".format(turningPoint,threshold))+\
    theme(figure_size=(16, 8))  # here you define the plot size+
    


p.save(filename = "{}/plot.png".format(outDir_COzones) )
logging.info("""    Files generated: 
    {}/plot.png """.format(outDir_COzones))

# %% End
logging.info("Finished")
