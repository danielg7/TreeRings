# Copyright notice: ----
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
# Script Intent: ---
# This model simulates the growth of three common savanna tree species
#
# Completeness: Incomplete
#
# Inputs: ----
# Add fuel moisture content.
#
#
# Outputs: ----
# 
# 
#
# TODO:  ----
#        
#        
#
#        
# 
# Load required packages ----

library(ggplot2)
library(ggthemes)

myTheme <- theme_tufte() +
  theme(
    text = element_text(family="Arial",size=18),
    axis.line = element_line(size = .3)
  )

#
# Run ancillary analyses:
source('/Users/danielgodwin/Dropbox/Graduate School/Dissertation/FireEscapeModel/Scripts/TreeHeight.R', echo=TRUE)


#
# Functions ----
savannaTree <- function(ID,
                        Species,
                        Height,
                        PersonalClock,
                        basalDiameter,
                        GrowthIntercept,
                        GrowthA,
                        GrowthError,
                        AnnualPrecipEst,
                        InteractEst,
                        TopKilled = FALSE)
{
  if(!is.numeric(Height) || !is.numeric(basalDiameter) || !is.numeric(GrowthIntercept) || !is.numeric(GrowthA) || !is.numeric(GrowthError)|| !is.numeric(AnnualPrecipEst)||!is.numeric(InteractEst))
     {
    stop("Height, basal diameter, and growth rate must all be numeric.")
  }
  st = data.frame(ID = ID,
            Species = Species,
            Height = Height,
            PersonalClock = PersonalClock,
            basalDiameter = basalDiameter,
            GrowthIntercept = GrowthIntercept,
            GrowthA = GrowthA,
            GrowthError = GrowthError,
            AnnualPrecipEst = AnnualPrecipEst,
            InteractEst = InteractEst,
            TopKilled = TopKilled)
  class(st) = "savannaTree"
  st
}

output_TreeMonitor <- function(savannaTree,TimeStamp,TimeSinceFire,monitorDF)
{
  scratchDF <- data.frame(ID = savannaTree$ID,
                          Species = savannaTree$Species,
                          Height = savannaTree$Height,
                          PersonalClock = savannaTree$PersonalClock,
                          basalDiameter = savannaTree$basalDiameter,
                          Time = TimeStamp,
                          RainfallRange = RainfallRange,
                          TopKilled = savannaTree$TopKilled,
                          TimeSinceFire = TimeSinceFire
                          )
  
  
  intermediateDF <- rbind(monitorDF,scratchDF)
  return(intermediateDF)
}

output_BurnMonitor <- function(Time,TimeSinceFire,RainfallRange,FRI,FuelLoad,Intensity,BurnMonitorDF)
{
  scratchDF <- data.frame(
    Time = Time,
                          TimeSinceFire = TimeSinceFire,
                          RainfallRange = RainfallRange,
                          FRI = FRI,
                          FuelLoad = FuelLoad,
                          Intensity = Intensity)
  intermediateDF <- rbind(BurnMonitorDF,scratchDF)
  return(intermediateDF)
}

probabilityOfTopKill <- function(ID,FireIntensity,TreeDiameter,TreeHeight,Method = "Diameter"){
  #if(!is.logical(UseDiameter)){stop("You must set whether to use diameter (default) or height for calculating topkill probability.")}
  
  if(Method == "Diameter"){
  topKillProb <- (4.5 -  3.4 * log(TreeDiameter) - 1.9)
  topKillProb <- exp(topKillProb)/(1+exp(topKillProb))}
  if(Method == "Intensity"){
    topKillProb <- 4.14e-05 * FireIntensity + 4.43e-01 }
  if(Method == "HeightAndIntensity"){
    topKillProb <- exp(4.3 - 5.003*log(TreeHeight) + 0.004408*sqrt(FireIntensity)) / (1 + exp(4.3 - 5.003*log(TreeHeight) + 0.004408*sqrt(FireIntensity)))
    #TopKillProb taken from Higgins et al. 2000 
  }
  
  topKillTest <- rbinom(1,1,topKillProb)
  
  returnValue = 1 #Default - no topkill
  
  if(topKillTest == 1)
  {
      
      returnValue = 0 #Return value indicates topkill
  }
  
  return(returnValue)
  }

    

  

#
# Initial Conditions ----
#

timeLimit <- 20
treeNumbers <- 1000
TimeSinceFire <- 0
PlotCounter <- 1
ResproutTrees <- TRUE
StochasticRainfall <- TRUE
RandomGrowth <- TRUE

RainfallRange <- c(250,550,750)
FRIRange <- c(1.1,2,3,4,5,6)

Combinations <- expand.grid(RainfallRange = RainfallRange,FRIRange = FRIRange)
# Empty dataframes for monitoring data----
monitorDF <- data.frame(ID = character(),
                        Species = character(),
                        Height = numeric(0),
                        PersonalClock = numeric(0),
                        basalDiameter = numeric(0),
                        Time = numeric(0),
                        RainfallRange = numeric(0),
                        FRI = numeric(0),
                        TopKilled = logical(0),
                        TimeSinceFire = numeric(0))

print(str(monitorDF))

BurnMonitor <- data.frame(
                          Time = numeric(0),
                          RainfallRange = numeric(0),
                          FRI = numeric(0),
                          Intensity = numeric(0))


metaMonitor <- monitorDF
# Main Script ----
for(m in 1:length(Combinations$FRIRange)){ # Script cycles through unique values in Fire Return Interval vector
  FireReturnInterval <- Combinations$FRIRange[m]
  FireFrequency <- 1/FireReturnInterval # FRI converted to fire frequency, to be used later in probability of a fire.

  

  if(StochasticRainfall == TRUE) {
    MAP <- rnorm(1,mean=Combinations$RainfallRange[m],sd=Combinations$RainfallRange[m] * 0.22294 + 107.76283)
    
    if(MAP <= 0){ 
      while(MAP <= 0) 
        MAP <- rnorm(1,mean=Combinations$RainfallRange[m],sd=Combinations$RainfallRange[m] * 0.22294 + 107.76283)}
  }
  if(StochasticRainfall == FALSE) MAP <- Combinations$RainfallRange[m]
  
for(k in seq(treeNumbers)) #Generate a cohort of trees
    {
    TESE <- savannaTree(ID = paste("TESE","-",FireReturnInterval,"-",MAP,"-","-",rnorm(1,100)), #Generates individuals in a cohort through a count of tree numbers, above.
                        Species = "TESE", #Assigns a species name
                        Height = .01, #Sets a base height (m)
                        PersonalClock = 0, 
                        basalDiameter = .1, #Sets a base basal diameter (cm)
                       # GrowthIntercept = -0.0439, #Sets the intercept
                      #  GrowthIntercept = coef(TESE_naiveGrowth)[1],
                        GrowthIntercept = 0,
                        GrowthA = coef(TESE_naiveGrowth)[1], #Sets the A of the model
                        GrowthError = summary(TESE_naiveGrowth)$coefficients[2,2] * sqrt(TESE_naiveGrowth$df.residual - 1),  #Sets the sd of the growth rate
                        AnnualPrecipEst = coef(TESE_naiveGrowth)[2],
                        InteractEst = 0,
                       )

    COAP <- savannaTree(ID = paste("COAP","-",FireReturnInterval,"-",MAP,"-","-",rnorm(1,100)),
                        Species = "COAP",
                        Height = .01,
                        PersonalClock = 0, 
                        basalDiameter = .1,
                        GrowthIntercept = 0,
                       # GrowthIntercept = coef(COAP_naiveGrowth)[1],
                        GrowthA = coef(COAP_naiveGrowth)[1],
                        GrowthError = summary(COAP_naiveGrowth)$coefficients[2,2] * sqrt(COAP_naiveGrowth$df.residual - 1),
                        AnnualPrecipEst = coef(COAP_naiveGrowth)[2],
                        InteractEst = 0)

    COMO <- savannaTree(ID = paste("COMO","-",FireReturnInterval,"-",MAP,"-","-",rnorm(1,100)),
                        Species = "COMO",
                        Height = .01,
                        PersonalClock = 0, 
                        basalDiameter = .1,
                        GrowthIntercept = 0,
                        #GrowthIntercept = coef(COMO_naiveGrowth)[1],
                        GrowthA = coef(COMO_naiveGrowth)[1],
                        GrowthError = summary(COMO_naiveGrowth)$coefficients[2,2] * sqrt(COMO_naiveGrowth$df.residual - 1),
                        AnnualPrecipEst = coef(COMO_naiveGrowth)[2],
                        InteractEst = 0)
    
    for(i in seq(1,timeLimit))
    { 
      if(ResproutTrees == TRUE){
        # Resprout if set to zero.
        if(COAP$basalDiameter <= 0){
          COAP$basalDiameter <- .1
          COAP$PersonalClock <- 0
        }
        
        if(TESE$basalDiameter <= 0)
        {
          TESE$basalDiameter <- .1
          TESE$PersonalClock <- 0
        }
        
        if(COMO$basalDiameter <= 0)
        {
          COMO$basalDiameter <- .1
          COMO$PersonalClock <- 0
        }
      }
      #Calculate height
      TESE$Height <- predict(TESE_nls,list(Diameter=TESE$basalDiameter))[1]
      COMO$Height <- predict(COMO_nls,list(Diameter=COMO$basalDiameter))[1]
      COAP$Height <- predict(COAP_nls,list(Diameter=COAP$basalDiameter))[1]
      
      # Figure out how much each tree would grow if it doesn't get top-killed
      
      if(RandomGrowth == TRUE){
     #   TESE_Growth_Step <- exp(TESE$GrowthIntercept + TESE$PersonalClock * rnorm(1,mean = TESE$GrowthA, sd = TESE$GrowthError) + TESE$AnnualPrecipEst * MAP + TESE$InteractEst * TESE$PersonalClock * MAP) 
      #  COAP_Growth_Step <- exp(COAP$GrowthIntercept + COAP$PersonalClock * rnorm(1,mean = COAP$GrowthA, sd = COAP$GrowthError) + COAP$AnnualPrecipEst * MAP + COAP$InteractEst * COAP$PersonalClock * MAP)
       # COMO_Growth_Step <- exp(COMO$GrowthIntercept + COMO$PersonalClock * rnorm(1,mean = COMO$GrowthA, sd = COMO$GrowthError) + COMO$AnnualPrecipEst * MAP + COMO$InteractEst * COMO$PersonalClock * MAP)
        
        COMO_Growth_Step <- predict(COMO_naiveGrowth_nls,list(Age=COMO$PersonalClock))[1]
        COAP_Growth_Step <- predict(COAP_naiveGrowth_nls,list(Age=COAP$PersonalClock))[1]
        TESE_Growth_Step <- predict(TESE_naiveGrowth_nls,list(Age=TESE$PersonalClock))[1]
        
      }
      if(RandomGrowth == FALSE){
        TESE_Growth_Step <- exp(TESE$GrowthIntercept + TESE$PersonalClock * TESE$GrowthA + TESE$AnnualPrecipEst * MAP + TESE$InteractEst * TESE$PersonalClock * MAP) 
        COAP_Growth_Step <- exp(COAP$GrowthIntercept + COAP$PersonalClock * COAP$GrowthA + COAP$AnnualPrecipEst * MAP + COAP$InteractEst * COAP$PersonalClock * MAP)
        COMO_Growth_Step <- exp(COMO$GrowthIntercept + COMO$PersonalClock * COMO$GrowthA + COMO$AnnualPrecipEst * MAP + COMO$InteractEst * COMO$PersonalClock * MAP)
      }
      #  TESE_Growth_Step <- exp(TESE$GrowthIntercept + TESE$PersonalClock * TESE$GrowthA + TESE$AnnualPrecipEst * MAP + TESE$InteractEst * TESE$PersonalClock * MAP) 
      #  COAP_Growth_Step <- exp(COAP$GrowthIntercept + COAP$PersonalClock * COAP$GrowthA + COAP$AnnualPrecipEst * MAP + COAP$InteractEst * COAP$PersonalClock * MAP)
      #  COMO_Growth_Step <- exp(COMO$GrowthIntercept + COMO$PersonalClock * COMO$GrowthA + COMO$AnnualPrecipEst * MAP + COMO$InteractEst * COMO$PersonalClock * MAP)
      
      
      # Did it burn?
      burnTest <- rbinom(n = 1,
                         size = 1,
                         prob = FireFrequency) #Probability is fire frequency, as defined by FRI in inits.
      if(burnTest == 0)
      {
        #Trees grow normally
        
        TESE$basalDiameter <- TESE$basalDiameter + TESE_Growth_Step
        COAP$basalDiameter <- TESE$basalDiameter + COAP_Growth_Step
        COMO$basalDiameter <- TESE$basalDiameter + COMO_Growth_Step
        
        TimeSinceFire <- TimeSinceFire + 1 #Add to a counter of time since fire if zero
        
        TESE$PersonalClock <- TESE$PersonalClock + 1
        COMO$PersonalClock <- COMO$PersonalClock + 1
        COAP$PersonalClock <- COAP$PersonalClock + 1
        
      }
      if(burnTest == 1)
      {
        FuelLoad <- 382.9 + 3.3 * MAP + 979.4 * TimeSinceFire - 0.001 * MAP^2 + 0.37*MAP*TimeSinceFire - 161.8*TimeSinceFire^2
        # Fuel loads calculated from Govender et al. 2006 Fig. 1
        #FuelLoad <- FuelLoad * .1 # Convert from kg / ha to g / m
        
        HeadFireOrBackfire <- rbinom(1,1,.7) #Randomly pick whether the plant will experience a headfire or backfire
        
        if(HeadFireOrBackfire == 1) #Values for headfire and backfire RoS are from Trollope 1978
          RateOfSpread <- rnorm(1,.15,.03)
        if(HeadFireOrBackfire == 0) 
          RateOfSpread <- rnorm(1,.02,.002)
        
        RelativeHumidity <- sample(seq(.042,.82,.01),1) # Range values taken from Trollope 2002
        FuelMoisture <- RelativeHumidity # Fuels assumed to take the value of RH immediately
        WindSpeed <- sample(seq(.3,6.7,.1),1) #(mean wind speed)
          
        #FireIntensity <- FuelLoad * 16.890 * RateOfSpread
        
        FireIntensity = 2729 + 0.8684*FuelLoad - 530*sqrt(FuelMoisture) - 0.1907*RelativeHumidity^2 - 5961/WindSpeed
        
        if(FireIntensity < 0){FireIntensity = 0}
        
        BurnMonitor <- output_BurnMonitor(Time = i,TimeSinceFire=TimeSinceFire,RainfallRange = MAP,FRI = FireReturnInterval,FuelLoad = FuelLoad,Intensity = FireIntensity,BurnMonitorDF = BurnMonitor)
        
        # Fire Intensity calculated as I = Hwr, where I is fire intensity (kW/m),
        # H is heat yield (kJ/g), w is the mass of fuel combused (g/m^2))
        # and r is the rate of spread (m/s)
        
        # Heat yield is derived from Govender et al. 2006
        # r is estimated from values seen in Govender et al. 2006
        TESE_Test <- probabilityOfTopKill(TESE$ID,FireIntensity,TESE$basalDiameter,TESE$Height,"HeightAndIntensity")
        if(TESE_Test == 0){
          TESE$TopKilled <- TRUE
          TESE$basalDiameter = 0
        }
        if(TESE_Test == 1){
          TESE$basalDiameter <- #TESE$basalDiameter + 
            TESE_Growth_Step
          TESE$TopKilled = FALSE
        }
        
        COAP_Test <- probabilityOfTopKill(COAP$ID,FireIntensity,COAP$basalDiameter,COAP$Height,"HeightAndIntensity")
        if(COAP_Test == 0){
          COAP$TopKilled <- TRUE
          COAP$basalDiameter = 0
        }
        if(COAP_Test == 1){
          COAP$basalDiameter <- #COAP$basalDiameter + 
            COAP_Growth_Step
          COAP$TopKilled <- FALSE}
        
        COMO_Test <- probabilityOfTopKill(COMO$ID,FireIntensity,COMO$basalDiameter,COMO$Height,"HeightAndIntensity")
        if(COMO_Test == 0){
          COMO$TopKilled <- TRUE
          COMO$basalDiameter <- 0}
        if(COMO_Test == 1){
          COMO$basalDiameter <- #COMO$basalDiameter + 
            COMO_Growth_Step
          COMO$TopKilled <- FALSE
        }
        
        TimeSinceFire <- 0
      }
      
      
      
     
      monitorDF <- output_TreeMonitor(savannaTree=TESE,i,TimeSinceFire,monitorDF)
      monitorDF <- output_TreeMonitor(savannaTree=COAP,i,TimeSinceFire,monitorDF)
      monitorDF <- output_TreeMonitor(savannaTree=COMO,i,TimeSinceFire,monitorDF)
    }

    
}    




monitorDF$RainfallRange <- Combinations$RainfallRange[m]
monitorDF$FRI <- FireReturnInterval
metaMonitor <- rbind(metaMonitor,monitorDF)

monitorDF <- data.frame(ID = character(),
                        Species = character(),
                        Height = numeric(0),
                        PersonalClock = numeric(0),
                        basalDiameter = numeric(0),
                        Time = numeric(0),
                        RainfallRange = numeric(0),
                        FRI = numeric(0))
}




metaMonitor$RainfallRange <- factor(metaMonitor$RainfallRange)
metaMonitor$FRI <- factor(metaMonitor$FRI)

metaMonitor <- metaMonitor[!duplicated(metaMonitor),]

# Plots ----


demoPlot <- ggplot(metaMonitor, aes(x = Time,
                                    y = Height,
                                    color = FRI,
                                    group = ID))
demoPlot+
  myTheme+
  scale_color_colorblind()+
  #geom_point(alpha=.5)+
  ylim(0,15)+
  geom_line(alpha=.25)+
  xlab("Time (Years...roughly)")+
  ylab("Height (m)")+
  facet_wrap(RainfallRange~Species,ncol=3)


metaMonitor_avg <- aggregate(data=metaMonitor,Height ~ Species + FRI + Time + RainfallRange,mean)

demoAvgPlot <- ggplot(metaMonitor_avg, aes(x = Time,
                                    y = Height,
                                    color = FRI))
demoAvgPlot+
  #ggtitle(paste("FRI ",m))+
  myTheme+
  scale_color_colorblind()+
  #geom_point(alpha=.5)+
  ylim(0,10)+
  geom_line(size=1.2,alpha=1)+
  xlab("Time (Years...roughly)")+
  ylab("Height")+
  facet_wrap(RainfallRange~Species,ncol=3)
# 
# 
# TopKillPlot <- ggplot(metaMonitor, aes(x = basalDiameter,
#                                        y = TopKillProb,
#                                        factor=RainfallRange,
#                                        color = Species, group=Species))
# TopKillPlot+
#   myTheme+
#   geom_point()+
#   geom_line()+
#   scale_color_discrete()+
#   xlab("Top Kill Probability")+
#   xlab("Basal Diameter (cm)")+
#   facet_wrap(~RainfallRange)



demoPlot_bwp <- ggplot(metaMonitor, aes(y = Height, x=FRI,factor=Species,fill=Species))
demoPlot_bwp+
  coord_flip()+
  #ggtitle(paste("FRI ",m))+
  myTheme+
  scale_fill_colorblind()+
  #geom_point(alpha=.5)+
  geom_boxplot()+
 # xlab("Time (Years...roughly)")+
#  ylab("Basal Diameter (cm)")+
  facet_wrap(~RainfallRange)

BurnMonitor$DidItBurn <- TRUE
BurnCheck <- merge(metaMonitor,BurnMonitor,by=c("Time","FRI","RainfallRange"),all.x=TRUE)
BurnCheck$DidItBurn <- replace(BurnCheck$DidItBurn,is.na(BurnCheck$DidItBurn),FALSE)

BurnCheck$DidItBurnGraph <- replace(BurnCheck$DidItBurn,BurnCheck$DidItBurn == TRUE,100)


BurnCheck$FRI <- as.numeric(BurnCheck$FRI)
BurnCheck$RainfallRange <- as.numeric(as.character(BurnCheck$RainfallRange))

BurnCheck <- subset(BurnCheck,Height <= 20)

demoPlotMulti <- ggplot(data = subset(BurnCheck,Species == "TESE"), aes(x = Time,
                                    y = Height,
                                    color = as.factor(FRI),
                                    group = ID))
demoPlotMulti+
  myTheme+
  scale_color_colorblind()+
  #geom_point(alpha=.5)+
  ylim(0,10)+
 #geom_point(data = subset(BurnCheck,!is.na(DidItBurnGraph) & Species=="TESE"),
  #          aes(x = Time, y = DidItBurnGraph),color="red",shape="2")+
  geom_line(alpha=.25)+
  xlab("Time (Years)")+
  ylab("Height (m)")+
  facet_wrap(FRI~RainfallRange,ncol=3)

BurnCheck$TopKilled <- replace(BurnCheck$TopKilled,BurnCheck$TopKilled==TRUE,1)
BurnCheck$TopKilled <- replace(BurnCheck$TopKilled,BurnCheck$TopKilled==FALSE,0)

BurnCounter <- ddply(BurnCheck,.(Species,RainfallRange,FRI,Height),summarise,
      CountKilled = sum(TopKilled == "1"),
      CountSurvive = sum(TopKilled=="0"))

BurnCounter$ProbTop <- BurnCounter$CountKilled / (BurnCounter$CountKilled + BurnCounter$CountSurvive)






testGLM <- glm(CountKilled ~ RainfallRange * Height * Species, data = BurnCounter, family=poisson)
testGLM_survive <- glm(CountSurvive ~ RainfallRange * Height * Species, data = BurnCounter, family=poisson)

#logitGLM <- glm(cbind(CountSurvive,CountKilled) ~ RainfallRange * Height * Species, data=BurnCounter,family=binomial(logit))

summary(testGLM)
summary(testGLM_survive)

test_wire <- expand.grid(Height = seq(0,10,2),
                       RainfallRange=seq(200,1000,100),
                       Species=c("COMO","TESE","COAP"))

 test_wire$TopKilled <- predict(testGLM,test_wire)
 test_wire$TopKilled_Survived <- predict(testGLM_survive,test_wire)
 test_wire$TopKilled <- exp(test_wire$TopKilled)
 test_wire$TopKilled_Survived <- exp(test_wire$TopKilled_Survived)
 
 test_wire$ProbTop <- test_wire$TopKilled / (test_wire$TopKilled+test_wire$TopKilled_Survived)
 
#test_wire$ProbTopLogit <- predict(logitGLM,test_wire)


newcols <- colorRampPalette(c("grey90", "grey10"))


wireframe(ProbTop ~ RainfallRange * Height | Species, data = test_wire,
          screen = list(z = 150, y = 0, x = -65),
          pretty = TRUE,
          drape = TRUE,
          scales = list(arrows=F,
                        cex=0.5
                        ),
          at = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),
          zlab = list("Probability of Top Kill",rot=90),
          xlab = list("\nMean Annual\nPrecipitation (mm)",rot=-20),
          ylab = list("Height (m)",rot=45),
          col.regions=newcols(100)
          )


write.csv(BurnCheck,file="Output/BurnCheck_1000.csv")
write.csv(BurnMonitor,file="Output/BurnMonitor_1000.csv")

