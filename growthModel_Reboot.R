# Copyright notice:
# This script is provided with a Creative Commons - Attribution license, as defined on:
# http://creativecommons.org/licenses/by/3.0/us/
#
#
# Author Contact:
# Daniel Godwin
# danielg7@gmail.com
# Savanna Ecology Lab
# Division of Biological Sciences
# University of Missouri - Columbia
#
# Script Intent:
# This model simulates the growth of three common savanna tree species
#
# Completeness: Incomplete
#
# Inputs:
#
#
#
# Outputs: 
# 
# 
#
# TODO:  Make growth a generalized function.
#        Make fire freq a fxn of rainfall
#        Make tree growth a function of rainfall
#        Make fire frequency, intensity a fxn of rainfall
# 
##################################
#
# Load required packages
#



library(ggplot2)
library(ggthemes)

myTheme <- theme_tufte() +
  theme(
    text = element_text(family="Arial",size=18),
    axis.line = element_line(size = .3)
  )

savannaTree <- function(ID,Species,Height,basalDiameter,GrowthIntercept,GrowthPolyA,GrowthPolyB,GrowthPolyC)
{
  if(!is.numeric(Height) || !is.numeric(basalDiameter) || !is.numeric(GrowthIntercept) || !is.numeric(GrowthPolyA) || !is.numeric(GrowthPolyB))
     stop("Height, basal diameter, and growth rate must all be numeric.")
  st = data.frame(ID = paste(Species,"-",ID,sep=""),
            Species = Species,
            Height = Height,
            basalDiameter = basalDiameter,
            GrowthIntercept = GrowthIntercept,
            GrowthPolyA = GrowthPolyA,
            GrowthPolyB = GrowthPolyB,
            GrowthPolyC = GrowthPolyC)
  class(st) = "savannaTree"
  st
}

output_TreeMonitor <- function(savannaTree,TimeStamp,monitorDF)
{
  scratchDF <- data.frame(ID = savannaTree$ID,
                          Species = savannaTree$Species,
                          Height = savannaTree$Height,
                          basalDiameter = savannaTree$basalDiameter,
                          Time = TimeStamp,
                          FRI = FireReturnInterval
                          )
  intermediateDF <- rbind(monitorDF,scratchDF)
  return(intermediateDF)
}

probabilityOfTopKill <- function(FireProb,FireIntensity,TreeDiameter){
  returnValue = 1
  
  if(TreeDiameter <= 0) {TreeDiameter == .1}
  
  if(rbinom(1,1,prob=FireProb) == 0)
    
  {
    topKillProb <- (4.5 -  3.4 * log(TreeDiameter) - 1.9)
    topKillProb <- exp(topKillProb)/(1+exp(topKillProb))
    
    
    topKillTest <- rbinom(1,1,topKillProb)
    
    if(topKillTest == 1)
    {
      returnValue = 0
    }
  }
  return(returnValue)
}



timeLimit <- 10 
treeNumbers <- 20
FireReturnInterval <- 2 #Annual Fire Return

FireFrequency <- 1 - 1/FireReturnInterval


monitorDF <- data.frame(ID = character(),
                        Species = character(),
                        Height = numeric(0),
                        basalDiameter = numeric(0),
                        Time = numeric(0),
                        FRI = numeric(0))

metaMonitor <- monitorDF

FireRange <- c(1,2,5,10)


for(p in unique(FireRange))
{
  FireReturnInterval <- k
  
for(k in seq(treeNumbers))
    {
    TESE <- savannaTree(k,"TESE",.1,.1,4.02415,39.50886,0.99700,1.80247)
    COAP <- savannaTree(k,"COAP",.1,.1,2.778,18.8894,-1.4695,0)
    COMO <- savannaTree(k,"COMO",.1,.1,2.42293,20.85294,5.61542,1.41780)
    
    monitorDF <- output_TreeMonitor(savannaTree=TESE,0,monitorDF)
    monitorDF <- output_TreeMonitor(savannaTree=COAP,0,monitorDF)
    monitorDF <- output_TreeMonitor(savannaTree=COMO,0,monitorDF)
    
    for(i in seq(1,timeLimit))
    { 
  
      TESE_Growth_Step <- TESE$GrowthIntercept + TESE$GrowthPolyA*i + TESE$GrowthPolyB*i^2 + TESE$GrowthPolyC*i^3
      COAP_Growth_Step <- COAP$GrowthIntercept + COAP$GrowthPolyA*i + COAP$GrowthPolyB*i^2 + COAP$GrowthPolyC*i^3
      COMO_Growth_Step <- COMO$GrowthIntercept + COMO$GrowthPolyA*i + COMO$GrowthPolyB*i^2 + COMO$GrowthPolyC*i^3
      
      TESE$basalDiameter <- TESE_Growth_Step * probabilityOfTopKill(FireFrequency,5,TESE$basalDiameter)
  
      if(TESE$basalDiameter == 0)
      TESE$basalDiameter <- .1

      COAP$basalDiameter <- COAP_Growth_Step * probabilityOfTopKill(FireFrequency,5,COAP$basalDiameter)
  
      if(COAP$basalDiameter == 0)
      COAP$basalDiameter <- .1
  
      COMO$basalDiameter <- COMO_Growth_Step * probabilityOfTopKill(FireFrequency,5,COMO$basalDiameter)
      
      if(COMO$basalDiameter == 0)
      COMO$basalDiameter <- .1
  
      monitorDF <- output_TreeMonitor(savannaTree=TESE,i,monitorDF)
      monitorDF <- output_TreeMonitor(savannaTree=COAP,i,monitorDF)
      monitorDF <- output_TreeMonitor(savannaTree=COMO,i,monitorDF)
  }
}

monitorDF$FRI <- p
metaMonitor <- rbind(metaMonitor,monitorDF)

}
##
#Plot the data
##

qplot(data = monitorDF,x = Time, y = basalDiameter, color = Species, geom = "point", facets = ~Species)

monitorDF$TopKillProb <- 4.5 -  3.4 * log(monitorDF$basalDiameter) - 1.9

monitorDF$TopKillProb <- exp(monitorDF$TopKillProb)/(1+exp(monitorDF$TopKillProb))

metaMonitor$TopKillProb <- 4.5 -  3.4 * log(metaMonitor$basalDiameter) - 1.9

metaMonitor$TopKillProb <- exp(metaMonitor$TopKillProb)/(1+exp(metaMonitor$TopKillProb))

metaMonitor$FRI <- factor(metaMonitor$FRI)

demoPlot <- ggplot(metaMonitor, aes(x = Time, y = basalDiameter, color = Species, linetype=ID))
demoPlot+
  myTheme+
  geom_point()+
  geom_line()+
  xlab("Time (Years...roughly)")+
  ylab("Basal Diameter (cm)")+
  facet_wrap(FRI~Species, scales="free_y",nrow=length(levels(metaMonitor$FRI)))



TopKillPlot <- ggplot(metaMonitor, aes(x = basalDiameter, y = TopKillProb, color = FRI, group=ID))
TopKillPlot+
  myTheme+
  geom_point()+
  geom_line()+
  scale_color_discrete()+
  xlab("Top Kill Probability")+
  xlab("Basal Diameter (cm)")+
  facet_wrap(FRI~Species, scales = "free_x",nrow=length(levels(metaMonitor$FRI)))


SizeClassPlotDF <- subset(metaMonitor,Time = max(Time))

SizeClassPlot <- ggplot(SizeClassPlotDF, aes(x = basalDiameter))
SizeClassPlot+
  myTheme+
  geom_histogram(binwidth = 20)+
  xlab("Basal Diameter (cm)")+
  facet_wrap(FRI~Species, nrow=length(levels(SizeClassPlotDF$FRI)))