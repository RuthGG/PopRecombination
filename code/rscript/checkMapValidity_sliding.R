#!/usr/bin/env Rscript
# Ruth Gómez Graciani
# 11 05 2022

###############################################################################
# Description:                                                                
### Map vs. map
# I compare ped/pop ratios (average, maximum and minimum) with expectations to make a p-value, and check how extreme are ratios compared to expected for each inversion.
# Sliding windows
# Inversions containing >80% of a window map
# 
# ### In vs. out
# I compare (average, maximum and minimum) in/out ratios with expectations to make p-value, and check how the ratio p-value is affected by inversion frequency on each map. 
# Sliding windows and raw data
# option 1: windows at breakpoints can be deleted, also 
# option 2: unless they overlap with >80% of either option, in which case they will be assigned to only in or only out

###############################################################################

# LOAD ARGUMENTS 
# print("### Load arguments")
# =========================================================================== #
args = commandArgs(trailingOnly=TRUE)

# # Test if there is at least one argument: if not, return an error
if (length(args)<7) {
  stop("Check inputs in script", call.=FALSE)
}

# # Example
# args[1]<-"analysis/defTest_LocationPatterns/04_MapValidityTest/inversionList/HsInv0371.bed" #Bedfile with invs to analyze
# args[2]<- 200000 #Confidence interval in bp
# args[3]<-1000 # Number of replicas/samples
# args[4]<-1
# args[5]<-"analysis/defTest_LocationPatterns/04_MapValidityTest/recombinationMap_sliding/fillWindowMeasurements/Spencemap_sliding/windowData.txt" # POP map in sliding windows made without repeat regions
# args[6]<-"analysis/defTest_LocationPatterns/04_MapValidityTest/recombinationMap_sliding/fillWindowMeasurements/Bherermap_sliding/windowData.txt" # PED map in sliding windows made without repeat regions
# args[7]<-outfile

# Process args/upload files
invs<-read.table(args[1], stringsAsFactors = F)
# just to make sure
invs<-invs[1,]
colnames(invs)<-c("Chromosome", "Start", "End", "inv")
CIsize<-as.numeric(args[2])
invs$CIStart<-invs$Start-CIsize
invs$CIEnd<-invs$End+CIsize
nSamples<-args[3]
requiredOverlap<-as.numeric(args[4])
genomePedname<-args[5]
genomePopname<-args[6]
outFile<-args[7]

# GENERAL FUNCTIONS

# A function that returns a subset with windows that have >=x% of overlap in a coordinate segment
# #             |====================|
# #   |-----------||----*---||---*-----|
# #   |------------------------------------|
takeInside<-function(what, from, overlap){ # Both in bed format 
  what[]<-lapply(what, as.character)
  what[,c(2,3)]<-lapply(what[,c(2,3)], as.numeric)
  from[]<-lapply(from, as.character)
  from[,c(2,3)]<-lapply(from[,c(2,3)], as.numeric)
  
  colnames(what)<-c("Chromosome", "Start", "End")
  colnames(from)<-c("Chromosome", "Start", "End")
  
  # Take everything that overlaps
  invPortion<-from[from$Chromosome == what$Chromosome & from$End >= what$Start & from$Start <= what$End,]
  
  # Make all starts and ends fit inside
  invPortion$Start_inside<-invPortion$Start
  invPortion$End_inside<-invPortion$End
  invPortion[invPortion$Start < what$Start, "Start_inside"]<- what$Start
  invPortion[invPortion$End > what$End, "End_inside"]<- what$End
  invPortion$size<-invPortion$End-invPortion$Start
  
  # Which is the window portion that is inside the inversion?
  invPortion$size_inside<-invPortion$End_inside-invPortion$Start_inside
  invPortion$overlap<-invPortion$size_inside/invPortion$size
  
  # Return a bedfile
  return(invPortion[invPortion$overlap >= overlap, c("Chromosome", "Start", "End")])
}


# INVERSION ASSESSMENT
# Load genomes
genomePed<-read.table(genomePedname, header = T, sep = "\t", stringsAsFactors = F)
genomePop<-read.table(genomePopname, header=T, sep ="\t", stringsAsFactors = F)

genomePed$ChrRegion<-genomePop$ChrRegion<-NULL
# We count how many windows fit inside this inv inv
# we want to make the comparison with the windows that have at least a 80% portion inside the inversion
invPortion<-takeInside(what =  invs[1,] , from =  genomePed, overlap = requiredOverlap)

winsInside<-nrow(invPortion)

if (winsInside){
  
  # MAKE SAMPLES
  # sampling function
  sampleGenome<-function(invInfo, genoInfo,  nSamples){
    ####################################################
    # invInfo is the inversion to simulate
    # genoInfo<-a recombination map, just to have the workspace available
    # nSamples<-number of samples we want to take
    #Make sure there is only 1 inv
    invInfo<-invInfo[1,]
    ####################################################
    
    # GENOME FILTERIGNG (remove inversion from sampleable space)
    # Tag genome
    rownames(genoInfo)<-paste(genoInfo$Chromosome, genoInfo$Start, genoInfo$End, sep = "_")
    genoInfo$ChrRegion<-"use"
    genoInfo$Chromosome<-as.character(genoInfo$Chromosome)
    genoInfo<-genoInfo[order(genoInfo$Chromosome, genoInfo$Start),]
    
    # Replace "use" with NA to avoid windows that  overlap with inv region
    # Any window overlapping with inv region and its confidence intervals, so not strict (those with *)
    #             |=================|   
    #   |-----*------||----*---||---*-----|
    #   |-------------------*-----------------|
    genoInfo[genoInfo$Chromosome == invInfo$Chromosome & genoInfo$End >= invInfo$CIStart & genoInfo$Start <= invInfo$CIEnd, "ChrRegion"]<-NA
    
    # TAKE SAMPLE COORDS
    # Now, a list of all the sampleable sites
    toSam<-genoInfo[which(genoInfo$ChrRegion == "use" ),]
    
    # We exclude coords at the end of usable fragments to avoid overhang
    toSam<-toSam[order(toSam$Chromosome, toSam$Start),]# sort
    toSam$previous<-toSam$Start - c(toSam$Start[2:length(toSam$Start)],0) # mark previous coord
    toSam$previousChr<-toSam$Chromosome == c(toSam$Chromosome[2:length(toSam$Chromosome)],0) # mark previous chrom
    # rownames(toSam)<-c(1:nrow(toSam)) # give rownames
    
    # Select unsampleable rownames
    slidingSize<-median(toSam$previous) # The most repeated number will be the sliding size
    
    nostarters<-toSam[which(toSam$previous != slidingSize | toSam$previousChr == FALSE),] #the chromosome is a double check...)
    
    nostarters_num<-which(rownames(toSam) %in% rownames(nostarters))
    nonums<-c()
    for(n in nostarters_num){
      nonums<-c(nonums, c((n-winsInside):n))
    }
    
    #Exclude coords again
    toSam[nonums, "ChrRegion"]<-NA
    toSam<-toSam[which(toSam$ChrRegion == "use" ),]
    
    # Make random samples
    samples<-dplyr::slice_sample(toSam, n = as.numeric(nSamples), replace = F)
    
    # Generate a table with all the info
    # samples<-toSam[starters, c("Chromosome", "Start")]
    samples$End<-samples$Start + (invInfo$End - invInfo$Start)
    samples$CIStart<-samples$Start - CIsize
    samples$CIEnd<-samples$End + CIsize
    
    #Add names
    samples$inv<-paste0(invInfo$inv, "_", c(1:nrow(samples)))
    return( samples[,c("Chromosome", "CIStart", "Start", "End", "CIEnd", "inv")])
    
    
  }
  # TAKE OVERLAPPING WINDOWS FROM SAMPLES
  fillSamples<-function(coords, genoInfo){
    
    # Unpack
    genoInfo$Chromosome<-as.character(genoInfo$Chromosome)
    
    # Find sliding windows for each section. We will use those windows overlapping with >80%, so 
    # breakpoints will have somewhat less coverage but that way we make sure that all included windows are representative
    
    # Wides
    winWide<-do.call("rbind", apply(coords,1,  function(x){
      invCoords<-data.frame(matrix(nrow = 1,c(as.character(x["Chromosome"]), as.numeric(x["CIStart"]), as.numeric(x["CIEnd"]))))
      inside<-takeInside(what = invCoords,from =  genoInfo, overlap = requiredOverlap)
      if(nrow(inside)>0){
        inside$inv<-x["inv"]
      }
      return(inside)
    }))
    winWide$name<-"out"
    # Insides
    winInside<-do.call("rbind", apply(coords,1,  function(x){
      invCoords<-data.frame(matrix(nrow = 1,c(as.character(x["Chromosome"]), as.numeric(x["Start"]), as.numeric(x["End"]))))
      inside<-takeInside(what = invCoords,from =  genoInfo, overlap = requiredOverlap)
      if(nrow(inside)>0){
        inside$inv<-x["inv"]
      }
      return(inside)
    }))
    winWide[rownames(winInside), "name"]<-"in"
    
    # Finish
    return(winWide)
    
  }
  
  samplesCoord<-sampleGenome(invs, genomePed, nSamples)
  samplesCoord$status<-"expected"
  invs$status<-"observed"
  samplesCoord<-rbind(samplesCoord, invs)
  expanded_samples<-fillSamples(samplesCoord, genomePed)
  
  # Add rates for both maps
  expanded_samples<-merge(expanded_samples, samplesCoord[,c("inv", "status")], by = "inv", all.x = T)
  samples_rated<-merge(expanded_samples, genomePed, by = c("Chromosome", "Start", "End"), all.x = T)
  samples_rated<-merge(samples_rated, genomePop,  by = c("Chromosome", "Start", "End"), all.x = T, suffixes = c(".ped", ".pop"))
  
  # MAP vs MAP
  # Make pedpop ratios
  samples_rated$ped.pop.ratio<-samples_rated$WAvgRate.ped/samples_rated$WAvgRate.pop
  
  names<-c("inv", "status")
  formula<-formula(paste0("ped.pop.ratio ~ ", paste(names, collapse = "+")))
  minmax<-  do.call(data.frame, aggregate( formula , samples_rated[samples_rated$name == "in",], FUN=function(x) c(min(x), max(x))))
  colnames(minmax)<-c(names, paste0( c("min", "max"), "PedPop_inside"))
  
  # I compare ped/pop ratios (maximum and minimum) with expectations to make a p-value.
  # Make a null distribution of expected min and max difference between maps in a sample this size
  nullRatioMin<-ecdf(minmax[minmax$status == "expected", "minPedPop_inside"]) # MIN
  nullRatioMax<-ecdf(minmax[minmax$status == "expected", "maxPedPop_inside"])# MAX
  
  # Compare observed to null distribution
  observedRatio<-minmax[minmax$status == "observed", ]
  observedRatio$minPedPop_inside_pval<-nullRatioMin(observedRatio$minPedPop_inside) # MIN
  observedRatio$maxPedPop_inside_pval<-nullRatioMax(observedRatio$maxPedPop_inside) # MAX
  
  # INSIDE vs OUTSIDE
  # I compare (average, maximum and minimum) in/out ratios with expectations to make p-value
  
  # Comparing function, gives MAX/MIN/AVG RATE INSIDE VS MAX/MIN/AVG RATE OUTSIDE 
  # Find sliding windows for samples (=inside + permutations) (InvsOutMeasurements )
  # Calculate min, max and avg recrate for outside & inside (InvsOutMeasurements)
  # For each measurement make the ratio between inside and outside (InvsOutMeasurements, column Fold_change)
  
  # Aggregate measurements (Calculate min, max and avg recrate for outside & inside )
  names<-c("name", "inv", "status")
  samples_rated_melt<-reshape2::melt(samples_rated, id.vars = c("Chromosome", "Start", "End", names))
  colnames(samples_rated_melt)[colnames(samples_rated_melt)== "variable"]<-"map"
  names<-c(names, "map")
  formula<-formula(paste0("value ~ ", paste(names, collapse = "+")))
  
  winMeasurements<-do.call(data.frame, (aggregate(formula, samples_rated_melt[samples_rated_melt$map != "ped.pop.ratio",], FUN = function(x) c(mean(x), min(x), max(x)))))
  colnames(winMeasurements)<-c(names, paste0("WAvgRate", c("", "_min", "_max")))
  
  # Make fold changes between measurements (in/out ratio)
  wMl<-reshape2::melt(winMeasurements, id.vars = names )
  wMw<-reshape(wMl, direction = "wide", idvar = c(names[names != "name"], "variable") , timevar = "name")
  wMw$Fold_change_inout<-wMw$value.in / wMw$value.out
  
  # Compare observed to null distribution for each map
  observedInout<-wMw[wMw$status == "observed",]
  expectedInout<-wMw[wMw$status == "expected",]
  
  observedInout$Fold_change_inout_pval<-apply(observedInout, 1, function(x){
    refTable<-expectedInout[expectedInout$map == as.character(x["map"]) & expectedInout$variable  ==  as.character(x["variable"]), "Fold_change_inout"]
    reference<-ecdf(refTable)
    reference(as.numeric(x["Fold_change_inout"]))
  })
  
  #MAKE A TABLE OF BOTH TOGETEHR
  final<-merge(observedRatio ,observedInout )
  final<-merge(final, invs, all.x = T)
  write.table(final, outFile)
  
}else{ 
  print(paste0("Inversion ", invs$inv, " did not have at least ",requiredOverlap*100,"% overlap with analysis windows"))
}
