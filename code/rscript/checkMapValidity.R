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
if (length(args)<6) {
  stop("Check inputs in script", call.=FALSE)
}

# # # Example
# args[1]<-"analysis/defTest_LocationPatterns/04_MapValidityTest/inversionList/test.bed" #Bedfile with invs to analyze
# args[2]<- 2000 #Confidence interval in bp
# args[3]<-10 # Number of replicas/samples
# args[4]<-"analysis/defTest_LocationPatterns/04_MapValidityTest/recombinationMap_sliding/fillWindowMeasurements/Spencemap_sliding/windowData.txt" # POP map in sliding windows
# args[5]<-"analysis/defTest_LocationPatterns/04_MapValidityTest/recombinationMap_sliding/fillWindowMeasurements/Bherermap_sliding/windowData.txt" # PED map in sliding windows
# args[6]<-outfile



# Process args/upload files
invs<-read.table(args[1], stringsAsFactors = F)
colnames(invs)<-c("Chromosome", "Start", "End", "inv")
CIsize<-as.numeric(args[2])
invs$CIStart<-invs$Start-CIsize
invs$CIEnd<-invs$End+CIsize
nSamples<-args[3]
genomePedname<-args[4]
genomePopname<-args[5]
outFile<-args[6]


# Make functions
sampleGenome<-function(invInfo, genoInfo,  nSamples){
  ####################################################
  # toSample is the list of invs we want to simulate
  # invInfo is the list of all invs (we want to avoid overlaps with these)
  # genome is the dataset to sample from
  # repDens can be "same" to resample in the main repDens profile of the real inv or "all" to resample in all the genome
  # nSamples is sample size invs, genomePed, nSamples
  # invInfo <- invs
  # genoInfo<-genomePed
  # # # repDens<-"same"
  # nSamples<-10
  #Make sure there is only 1 inv
  invInfo<-invInfo[1,]
  ####################################################
  
  # GENOME FILTERIGNG
  
  # Tag genome
  rownames(genoInfo)<-paste(genoInfo$Chromosome, genoInfo$Start, genoInfo$End, sep = "_")
  genoInfo$ChrRegion<-"use"
  # Replace Repeat info with NA to avoid windows that  overlap with inv region
  genoInfo[genoInfo$Chromosome == invInfo$Chromosome & genoInfo$Start >= invInfo$CIStart & genoInfo$End <= invInfo$CIEnd, "ChrRegion"]<-NA
  
  # We count how many invs do we need for this inv
  sampleSize<-sum(is.na(genoInfo$ChrRegion))
  
  # Now, for each used repeat status OR for all the genome, we want to make a list of all the sampleable sites
  genoInfo<-genoInfo[!is.na(genoInfo$ChrRegion) ,]
  
  # We exclude coords at the end of usable fragments to avoid overhang
  toSam<-genoInfo
  toSam<-toSam[order(toSam$Chromosome, toSam$Start),]# sort
  toSam$previous<-toSam$Start - c(0, toSam$Start[1:length(toSam$Start)-1]) # mark previous coord
  rownames(toSam)<-c(1:nrow(toSam)) # give rownames
  
  # Select sampleable rownames
  slidingSize<-median(toSam$previous) # The most repeated number will be the sliding size
  nostarters<-c(as.numeric(rownames(toSam[toSam$previous != slidingSize,])), nrow(toSam)+1)-sampleSize+1
  nonums<-c()
  for(n in nostarters){
    nonums<-c(nonums, c(n:(n+sampleSize-2)))
  }
  
  # Make a list of options
  options<-as.numeric(rownames(toSam)[!(as.numeric(rownames(toSam)) %in% nonums)])
  
  # Make random samples
  starters<-sample(options, nSamples, replace = F)
  
  # Generate a table with all the info
  samples<-toSam[starters, c("Chromosome", "Start")]
  samples$End<-samples$Start + (invInfo$CIEnd - invInfo$CIStart)
  colnames(samples)<-c("Chromosome", "CIStart", "CIEnd")
  samples$Start<-samples$CIStart + CIsize
  samples$End<-samples$CIEnd -CIsize
  
  #Add names
  samples$inv<-paste0(invInfo$inv, "_", c(1:nrow(samples)))
  return(samples)
}
InvsOutMeasurements<-function(table, genoInfo){
  ##################
  # table<-samples
  # genoInfo<-genome
  ##################
  
  # Unpack
  newcols<-c("Chromosome", "Start", "End", "inv")
  left<-table[,c("Chromosome", "CIStart", "Start", "inv")]
  colnames(left)<-newcols
  left$name<-"out"
  right<-table[,c("Chromosome", "End", "CIEnd", "inv")]
  colnames(right)<-newcols
  right$name<-"out"
  inv<-table[,c("Chromosome" , "Start", "End", "inv")]
  colnames(inv)<-newcols
  inv$name<-"in"
  sections<-rbind(rbind(left, right), inv)
  rm(left, right, inv)
  
  # Find sliding windows for samples
  winSections<-do.call("rbind", apply(sections,1,  function(x){
    b<-genoInfo[genoInfo$Chromosome == as.character(x["Chromosome"]) & (genoInfo$Start >= as.numeric(x["Start"]) & genoInfo$End <= as.numeric(x["End"])),]
    if(nrow(b)>0){
      b$name<-x["name"]
      b$inv<-x["inv"]
    }
    return(b)
  }))
  
  # Aggregate measurements (Calculate min, max and avg recrate for outside & inside )
  names<-c("name", "inv")
  formula<-formula(paste0("WAvgRate ~ ", paste(names, collapse = "+")))
  
  avg<-aggregate(formula, winSections, "mean")
  min<-aggregate(formula, winSections, "min")
  max<-aggregate(formula, winSections, "max")
  
  winMeasurements<-merge(merge(max, min, by = names, suffixes = c("_max", "_min")), avg, by = names)
  
  # Make fold changes between measurements
  wMl<-reshape2::melt(winMeasurements, id.vars = c("name", "inv") )
  wMw<-reshape(wMl, direction = "wide", idvar = c("inv", "variable"), timevar = "name")
  wMw$Fold_change<-wMw$value.in / wMw$value.out
  return(wMw)
}
MapvsMapMeasurements<-function(table, genoInfoA, genoInfoB){
  
  #################
  # genoInfoA<-genomePed
  # genoInfoB<-genomePop
  # table<-samplesCoord
  # x<-table[1,]
  ################
  table$Chromosome<-as.character(table$Chromosome)
  table[,c("minFC", "maxFC")]<-    apply(table,1,  function(x){
    a<-genoInfoA[genoInfoA$Chromosome == as.character(x["Chromosome"]) & (genoInfoA$Start >= as.numeric(x["Start"]) & genoInfoA$End <= as.numeric(x["End"])),"WAvgRate"]
    b<-genoInfoB[genoInfoB$Chromosome == as.character(x["Chromosome"]) & (genoInfoB$Start >= as.numeric(x["Start"]) & genoInfoB$End <= as.numeric(x["End"])),"WAvgRate"]
    sectionpedpopratios<-a/b
    return(c(min(sectionpedpopratios), max(sectionpedpopratios)))
  })
  
  return(table)
  
}


# QUANTILE OF MIN AND MAX RATIO BETWEEN MAPS INSIDE INV

genomePed<-read.table(genomePedname, header = T)
genomePop<-read.table(genomePopname, header=T)
# Make permutations (sampleGenome)
samplesCoord<-sampleGenome(invs, genomePed, nSamples)
# Find sliding windows for samples (=inside + permutations)
# For each window calculate ped/pop ratio
# Find min and max for the ratio in each sample
samplesRatio<-MapvsMapMeasurements(samplesCoord, genomePed, genomePop)
observedRatio<-MapvsMapMeasurements(invs, genomePed, genomePop)
# Make a null distribution of expected min and max difference between maps in a sample this size
nullRatioMin<-ecdf(samplesRatio$minFC)
nullRatioMax<-ecdf(samplesRatio$maxFC)
# Compare observed to null distribution
observedRatio$minFC_pval<-nullRatioMin(observedRatio$minFC)
observedRatio$maxFC_pval<-nullRatioMax(observedRatio$maxFC)

#COMPARISON OF QUANTILES OF THE DIFFERENCE in various measurements BETWEEN INSIDE AND OUTSIDE

# Make permutations (sampleGenome) -> samplesCoord
# Find sliding windows for samples (=inside + permutations) (InvsOutMeasurements )
# Calculate min, max and avg recrate for outside & inside (InvsOutMeasurements)
# For each measurement make the ratio between inside and outside (InvsOutMeasurements, column Fold_change)
# Make a null distribution of the expected diff inside vs outside for each map
samplesInoutPed<-InvsOutMeasurements(samplesCoord, genomePed )
samplesInoutPop<-InvsOutMeasurements(samplesCoord, genomePop)
# Compare observed to null distribution for each map
observedInoutPed<-InvsOutMeasurements(invs, genomePed)
observedInoutPed$Fold_change_pval<-apply(observedInoutPed, 1, function(x){
  reference<-ecdf(samplesInoutPed[samplesInoutPed$variable ==  as.character(x["variable"]), "Fold_change"])
  reference(as.numeric(x["Fold_change"]))
})
observedInoutPop<-InvsOutMeasurements(invs, genomePop)
observedInoutPop$Fold_change_pval<-apply(observedInoutPop, 1, function(x){
  reference<-ecdf(samplesInoutPop[samplesInoutPop$variable ==  as.character(x["variable"]), "Fold_change"])
  reference(as.numeric(x["Fold_change"]))
})
# Compare quantiles between both maps
observedInout<-merge(observedInoutPed, observedInoutPop, by = c("inv", "variable"), suffixes = c("_ped", "_pop"))
observedInout$pvalRatio<-observedInout$Fold_change_pval_ped / observedInout$Fold_change_pval_pop
observedInout$pvalDiff<-observedInout$Fold_change_pval_ped - observedInout$Fold_change_pval_pop


#MAKE A TABLE OF BOTH TOGETEHR
final<-merge(observedRatio , observedInout[,c("inv", "variable", "Fold_change_ped", "Fold_change_pop", "Fold_change_pval_ped", "Fold_change_pval_pop", "pvalRatio", "pvalDiff")])
write.table(final, args[6])

