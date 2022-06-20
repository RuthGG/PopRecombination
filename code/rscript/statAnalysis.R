#!/usr/bin/env Rscript
# Ruth Gómez Graciani
# 11 05 2022

###############################################################################
# Description:                                                                
# Make statistical analysis for a directory      
###############################################################################

# LOAD ARGUMENTS 
print("### Load arguments")
# =========================================================================== #
args = commandArgs(trailingOnly=TRUE)

# # Test if there is at least one argument: if not, return an error
if (length(args)<4) {
  stop("One path option within divideChromosomes, one path option within fillWindowMeasurements, one output directory, one option name", call.=FALSE)
}

# # Example
# args[1]<-"analysis/20220511_LocationPatterns/divideChromosomes/avgBherer_fixedArms_5/" # Directory for divideChromosomes
# args[2]<-"analysis/20220511_LocationPatterns/fillWindowMeasurements/avgBherer_fixedArms_5/"  # Directory with fillWindowMeasurements
# args[3]<-"analysis/20220511_LocationPatterns" # Exit file
# args[4]<-"avgBherer_fixedArms_5/"

# LOAD PACKAGES
print("### Load packages")
# =========================================================================== #
if (!require("pacman")) install.packages("pacman")
divChroms<- args[1]
winMeasur<- args[2]
dir.create(paste0(args[3], "/statAnalysis/"))
outdir<- paste0(args[3], "/statAnalysis/",args[4])
dir.create(outdir)
name<-args[4]

### Workspace visualization
print("### Workspace visualization")
# =========================================================================== #

pacman::p_load(GenomicRanges, ggplot2, ggbio)
data(ideoCyto, package = "biovizBase")

visWorksp <- function(dataLimitsFile, name){
  dataLimits <- read.table(dataLimitsFile ,sep = "\t", header = T)
  
  # ARM REGIONS 
  arms.GR<-makeGRangesFromDataFrame(dataLimits[grep("cen", dataLimits$chromID, invert = T),], ignore.strand = T, start.field = "Start", end.field = "End" , seqnames.field = "Chromosome")
  
  # ORIGIN OF LIMITS
  viewWidth = 100
  dataLimits$Start_data_end <- dataLimits$Start + viewWidth
  dataLimits$End_data_start <- dataLimits$End - viewWidth
  
  limits.GR <-makeGRangesFromDataFrame(
    rbind(
      setNames(data.frame(dataLimits[ ,c("Chromosome", "Start","Start_data_end")]), c("chrom","Start","End")) ,
      setNames(data.frame(dataLimits[,c("Chromosome", "End_data_start", "End")]), c("chrom","Start", "End"))
    ),
    keep.extra.columns = TRUE)
  
  # PLOT 
  autoplot(ideoCyto$hg19, layout = "karyogram", cytobands = TRUE)+
    layout_karyogram(data = arms.GR, geom = "rect", ylim = c(0, 10), fill = "yellow", alpha = 0.3)+
    layout_karyogram(data = limits.GR, geom = "rect", ylim = c(0, 10), color = "orange",  alpha = 1, size = 1)+
    guides(fill = "none")+
    ggtitle(name)
}
visWindows <- function(divChroms, name){
  
  # PATHS
  windowsFile <- paste0(divChroms,"/windows.txt")
  densityFile <- paste0(divChroms,"/densities.txt")
  extremesFile <- paste0(divChroms,"/extremes.txt")
  centroFile <- paste0(divChroms,"/workspace.txt")
  
  # FILES
  windows <- read.table(windowsFile, header = T, sep = "\t")
    windows$Color <- "a"
    windows[as.numeric(rownames(windows)) %% 2 == 1,"Color"]<-"b"
  windows$Chromosome <- factor(windows$Chromosome, levels =  paste(rep("chr", 23), as.character(c(c(1:22),"X")), sep = ""))
    
  density <- read.table(densityFile, header = T, sep = "\t")
  density$Chromosome <- factor(density$Chromosome, levels = paste(rep("chr", 23), as.character(c(c(1:22),"X")), sep = ""))
  
  extremes <- read.table(extremesFile, header = T, sep = "\t")
  extremes$Chromosome <- factor(extremes$Chromosome, levels = paste(rep("chr", 23), as.character(c(c(1:22),"X")), sep = ""))
  
  centromeres <- read.table(centroFile, header = T, sep = "\t")
  centromeres$Chromosome <- factor(centromeres$Chromosome, levels =  paste(rep("chr", 23), as.character(c(c(1:22),"X")), sep = ""))
  armLimits <-centromeres[grep("cen", centromeres$chromID,invert = T),]
  centromeres<-centromeres[grep("cen", centromeres$chromID),]
  starts <-  armLimits[grep("p", armLimits$chromID), c("Start", "Chromosome")]
  ends<- armLimits[grep("q", armLimits$chromID), c("End", "Chromosome")]
  
  # PLOT
  ggplot()+
        geom_rect(data= centromeres,  aes(xmin = Start, xmax = End,  ymin = 0, ymax = Inf), fill = "blue4", alpha = 0.5)+
        geom_vline(data=starts, aes(xintercept = Start), color = "orange")+geom_vline(data=ends, aes(xintercept = End), color = "orange")+
        geom_rect(data=windows, aes(xmin = Start, xmax = End, fill = Color, ymin = 0, ymax = Inf), alpha = 0.3)+
        geom_line(data=density,aes(x = pos, y = val))+
        geom_point(data=extremes, aes(x = pos, y = val, color = Type))+
        facet_wrap("Chromosome", scales = "free")+
        scale_fill_manual(values=c("#737373", "#e1e5eb"), guide="none")+
        scale_color_manual(values = c("#bd2b43", "#2b63bd"), guide = "none")+
        ggtitle(name)
  
}

# visWorksp(  dataLimitsFile = "analysis/20220505_LocationPatterns/divideChromosomes/mapChromosomeLimits_avgBherer.txt", "Bherer map limits")
# visWindows(divChroms = "analysis/20220511_LocationPatterns/divideChromosomes/avgBherer_COzones_0_800000/", "Crossover zones for Bherer map")


# Detect analysis type and write figure

if( file.exists(paste0(divChroms,"/densities.txt")) ){
  figscale= 2
  png(paste0(outdir, "/cozones.png"), width=2480/figscale, height=1748/figscale, units = "px")
  print(visWindows(divChroms = divChroms, paste("Chrossover zones for", name)))
  dev.off()
}else{
  figscale = 3
  png(paste0(outdir, "/workspace.png"), width=2480/figscale, height=1748/figscale, units = "px")
  print(visWorksp(dataLimitsFile = paste0(divChroms,"/workspace.txt"), paste("Map limits (workspace) for", name)))
  dev.off()
}


### Variable visualization
print("### Variable visualization")
# =========================================================================== #

pacman::p_load(ggdist, ggplot2, gghalves, reshape2, patchwork)

visVars <- function(windowDatafile, name){
  # Read file
  windowData <- read.table(windowDatafile, sep = "\t", header = T)
  
  # Put chromosome factor in order
  windowData$Chromosome <- factor(windowData$Chromosome, levels = paste(rep("chr", 22), as.character(c(1:22)), sep = ""))
  windowData$chromID <- windowData$winID <- NULL
  
  # Melt data
  windowDataMelted <- melt(windowData, id.vars = c( "Chromosome", "Start", "End")) 
  
  # Group info
  windowDataMelted$dataGroup<- ifelse(windowDataMelted$variable %in% c("invCenters", "NHCenters", "NAHRCenters"), "Inversions", 
                                      ifelse(windowDataMelted$variable %in% c("allRepCounts", "intraRepCounts"), "Repeats",
                                             ifelse(windowDataMelted$variable %in% c("WAvgRate"), "Weighted average recRate",
                                                    ifelse(windowDataMelted$variable %in% c("maxRate"), "Maximum recRate",  "Window length"))))
  
  # Add log repeats
  windowDataMelted_replog <- windowDataMelted[windowDataMelted$dataGroup == "Repeats",]
  windowDataMelted_replog$value <- log10(windowDataMelted_replog$value)
  windowDataMelted_replog$dataGroup <- "Repeats_log10"
  windowDataMelted<- rbind(windowDataMelted, windowDataMelted_replog)
  
  # Make list of plots
  plot_list<-list()
  for (group in unique(windowDataMelted$dataGroup)) {
    
    plotTable <- windowDataMelted[(windowDataMelted$value != -Inf) & (windowDataMelted$dataGroup == group),]
    plot_list[[group]] <- ggplot(plotTable, aes(x = variable, y = value))+
      # Half violin
      ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.2, point_colour = NA) +
      # Boxplot 
      geom_boxplot(width = .1, outlier.shape = NA) +
      # Points
      gghalves::geom_half_point_panel(side = "l", range_scale = .6,  alpha = .5, aes(color = Chromosome))+
      scale_color_manual(values = c(rep("#3c7ae7",11),rep("#89b23e",11) ))+
      # Adjust coordinates
      coord_flip()+
      # coord_flip( xlim = c(1.3, NA))+
      # Adjust labels
      theme(axis.title.y = element_blank(), legend.position = "none")+
      # Title
      ggtitle(group)
  }
  
  # Plot list of plots
  wrap_plots(plot_list)+ plot_annotation( title = name)
}

# visVars("analysis/20220509_LocationPatterns/fillWindowMeasurements/windowData_0_800000_avgBherer.txt", "avgBherer")

# Detect analysis type and write figure
figscale = 2
  png(paste0(outdir, "/variableVisualization.png"), width=2480/figscale, height=1748/figscale, units = "px")
  print(visVars(windowDatafile =  paste0(winMeasur,"/windowData.txt"), paste("Variable distributions for", name)))
  dev.off()

### Correlations between variables
print("### Correlations between variables")
# =========================================================================== #

pacman::p_load(corrplot, PerformanceAnalytics)

figscale = 2

  # Read file
  windowDatafile <-  paste0(winMeasur,"/windowData.txt")
  windowData <- read.table(windowDatafile, sep = "\t", header = T)
  windowData$chromID <- windowData$winID <- NULL
  
  # Make new table
  corrTable<-windowData[,!(colnames(windowData) %in% c("Start", "End", "Chromosome"))]
  
  # Make correlation plots
  png(paste0(outdir, "/correlationVisualization_Pearson.png"), width=2480/figscale, height=1748/figscale, units = "px")
  chart.Correlation(corrTable, histogram=TRUE, pch=19, method = "pearson")
  mtext(paste0(name, ". Pearson correlation"), line = 3)
  dev.off()
  
  png(paste0(outdir, "/correlationVisualization_Spearman.png"), width=2480/figscale, height=1748/figscale, units = "px")
  chart.Correlation(corrTable, histogram=TRUE, pch=19, method = "spearman")
  mtext(paste0(name, ". Spearman correlation"), line = 3)
  dev.off()
  
## Multiple linear regression
print("### Multiple linear regression")
# =========================================================================== #

pacman::p_load(MASS, performance)

makeModel<-function(windowDatafile,var){
  # Read  file
  windowData <- read.table(windowDatafile, sep = "\t", header = T)
  
  # Make total model
  allvars<-lm(paste0(var, " ~   allRepCounts + WAvgRate +  Length.bp.  + maxRate"), data = windowData)
  outallvars<-capture.output(summary(allvars))
  
  # Check optimal model
  step <- stepAIC(allvars, direction="both")
  outstep<-capture.output(step$anova)
  
  # Read  file... again
  windowData <- read.table(windowDatafile, sep = "\t", header = T)
  
  # Make optimal model
  stepvars<-lm(formula(step), data = windowData)
  outstepvars<-capture.output(summary(stepvars))
  
  # Write all outputs
  info <- paste0(c(paste0("Checking model for variable: ",var), "#================================#","#---------------------------------#",
                   "All-in model", "#---------------------------------#",  outallvars,"#---------------------------------#",
                   "Step fitting", "#---------------------------------#", outstep, "#---------------------------------#",
                   "Final model", "#---------------------------------#", outstepvars),
                 collapse="\n")
  
  # Evaluate optimal model
  plot<- check_model(stepvars)
  
  #Return
  return(list(info = info,plot=plot))
}

# test<-makeModel("analysis/20220509_LocationPatterns/fillWindowMeasurements/windowData_0_800000_avgBherer.txt","NAHRCenters" )

# Setup
vars = c("invCenters", "NHCenters", "NAHRCenters")
figscale = 2

# Run
for(var in vars){
    dirsave <- paste0(outdir,"/model_",var)
    dir.create(dirsave)
    
    mod<- makeModel(paste0(winMeasur,"/windowData.txt") ,var )
    
    cat(mod$info, file = paste0(dirsave, "/stepModel.txt"), append = F)
    
    png( paste0(dirsave, "/evalModel.png"), width=2480/figscale, height=1748/figscale, units = "px")
    print(mod$plot)
    dev.off()
}


