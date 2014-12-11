timeLength <- 500

####
# Create dummy lookup table
###

TESE_Growth_lower <- data.frame(Species = "TESE",
                          basalDiameterIncrement = 1.05,
                          heightIncrement = 1.08,
                          minimumRainfall = 700,
                          maximumRainfall = 800)

TESE_Growth_upper <- data.frame(Species = "TESE",
                          basalDiameterIncrement = 1.07,
                          heightIncrement = 1.10,
                          minimumRainfall = 800,
                          maximumRainfall = 1000)


lookupTable <- rbind(TESE_Growth_upper,TESE_Growth_lower)

growthRate <- function(SpeciesType, basalDiameter, Height, annualRainfall, returnType)
{
  if(returnType == "diameterIncrement")
  {
    lookupTable_subset <- subset(lookupTable, lookupTable$Species == SpeciesType)
    lookupTable_subset <- subset(lookupTable_subset, annualRainfall >= lookupTable_subset$minimumRainfall & annualRainfall <= maximumRainfall )
    
    diameterIncrement <- basalDiameter * lookupTable_subset$basalDiameterIncrement
    
    return(diameterIncrement)
  }
  
if(returnType == "heightIncrement")
{
  lookupTable_subset <- subset(lookupTable, lookupTable$Species == SpeciesType)
  lookupTable_subset <- subset(lookupTable_subset, annualRainfall >= lookupTable_subset$minimumRainfall & annualRainfall <= maximumRainfall )
  
  heightIncrement <- Height * lookupTable_subset$heightIncrement

return(heightIncrement)
}
}

treeSurvival <- data.frame(Species = character(0),
                           basalDiameter = numeric(0),
                           Height = numeric(0),
                           Time = numeric(0))

for( i in seq(from = 1, to = timeLength))
{
  
  scratchSurvival <- data.frame(Species = "TESE",
                                basalDiameter = growthRate(SpeciesType = "TESE",
                                                           basalDiameter = treeSurvival[i,2],
                                                           Height = treeSurvival[i,3],
                                                          annualRainfall = 850,
                                                          returnType = "basalIncrement"),
                                Height = growthRate(SpeciesType = "TESE",
                                                          basalDiameter = treeSurvival[i,2],
                                                          Height = treeSurvival[i,3],
                                                           annualRainfall = 850,
                                                           returnType = "heightIncrement"),
                                Time = i) 
  # * FireEffect()
  
  treeSurvival <- rbind(treeSurvival,scratchSurvival)
}


# TreeGrowth * ProbabilityOfTopKill

