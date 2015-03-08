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
# This script combines known annual rainfall data from Kruger National Park rain gauges
# and uses it to estimate and fill in the gaps.
#
# Completeness: Incomplete
#
# Inputs: ----
# Kruger_Stations1980_2008.csv
# monthly_rainfall_2006_2010_paf-moo-sat-pre.txt
# WRF_rainfall_annual_totals.csv
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

library(reshape)
library(plyr)
library(ggplot2)
library(lattice)
library(corrplot)
library(dplyr)

setwd("/Users/danielgodwin/Dropbox/Graduate School/Dissertation/Chapter 3 - Growth/TreeRings/")

sitesWx <- read.csv("Weather/Kruger_Stations1980_2008.csv")

sitesWx_long <- melt(sitesWx,id.vars="Year",variable_name="Station")
names(sitesWx_long)[3] <- "AnnualPrecip"

newWx <- read.csv("Weather/monthly_rainfall_2006_2010_paf-moo-sat-pre.txt")
newWx_annual <- ddply(newWx,.(YEAR,STATION),summarise,AnnualPrecip = sum(SumOfMM))
names(newWx_annual) <- c("Year","Station","AnnualPrecip")

newWx_annual_subset <- subset(newWx_annual,Year >= 2008)

Kruger_Wx_Combined <- rbind(sitesWx_long,newWx_annual_subset)

lookupStation <- data.frame(Station = unique(Kruger_Wx_Combined$Station),Station_Long = c("Satara","BergEnDal","Byamiti","Houtboschrand","Kingfisherspruit","CrocodileBridge","Letaba","Mahlangeni","Malelane","Mooiplaas","Mopani","Nwanetsi","Olifants",
                                                                                          "Lower Sabie","Pafuri","PafuriWenela","Phalaborwa","Pretoriuskop","PundaMaria","Shangoni","Shingwedzi","Shimuwini","Sirheni","Skukuza","Stolznek","Talamati",
                                                                                          "Tshokwane","Vlakteplaas","Woodlands"))

Kruger_Wx_Combined <- merge(Kruger_Wx_Combined,lookupStation,by="Station")


WRF_Wx <- read.csv("Weather/WRF_rainfall_annual_totals.csv")
names(WRF_Wx)[2] <- "AnnualPrecip"

WRF_Wx$Station <- "WRF"
WRF_Wx$Station_Long <- "WitsRuralFacility"

Kruger_Wx_Combined <- rbind(Kruger_Wx_Combined,WRF_Wx)


# Do correlation table ----------------------------------------------------

corrTableWx_wide <- dcast(data = Kruger_Wx_Combined,formula = Year ~ Station_Long, max, value.var = "AnnualPrecip")
corrTableWx_wide$Year <- NULL
corrTableWx_wide <- do.call(data.frame,lapply(corrTableWx_wide, function(x) replace(x, is.infinite(x),NA)))

df2 <- cor(corrTableWx_wide, use = "na.or.complete")
corrplot(df2, method="shade",shade.col=NA, tl.col="black", tl.srt=45)

MakeJustWide <- function(DF = data.frame(), GoodStation, BadStation){
  df_intermediate <- subset(DF,Station_Long == BadStation | Station_Long == GoodStation)
  qplot(data = df_intermediate, x = Year, y = AnnualPrecip, color = Station_Long, geom = c("line","point"))+theme_bw()
  df_wide <- dcast(data = df_intermediate,formula = Year ~ Station_Long, mean, value.var = "AnnualPrecip")
  print(qplot(data = df_wide, x = df_wide[,2], y = df_wide[,3],xlab = GoodStation, ylab = BadStation))
  
  return(df_wide)
}


# BergEnDal Corrections ---------------------------------------------------

# Fill in Berg En Dal. Malelane is the closest station.
# Subset first
a <- MakeJustWide(Kruger_Wx_Combined,GoodStation = "Malelane", BadStation = "BergEnDal")

# Test correlation
cor.test(x = a$BergEnDal,y = a$Malelane)

# Fit linear model
Berg_By_Mal <- lm(data = a, BergEnDal ~ Malelane)
summary(Berg_By_Mal)

# Find the holes
Berg_and_Mal_Wx_holes <- subset(a,is.na(BergEnDal))

# Fill the holes
Berg_and_Mal_Wx_holes$BergEnDal <- predict(object = Berg_By_Mal,newdata = Berg_and_Mal_Wx_holes)

# Check to make sure it looks good
qplot(data = a, x = Malelane, y = BergEnDal) + geom_point(data = Berg_and_Mal_Wx_holes, colour="red", aes(x = Malelane, y = BergEnDal))

rm(a)

# Byamiti Corrections ---------------------------------------------------

# Fill in Byamiti. Malelane is the closest station.
# Subset first
a <- MakeJustWide(Kruger_Wx_Combined,GoodStation = "CrocodileBridge", BadStation = "Byamiti")

# Test correlation
cor.test(x = a$CrocodileBridge,y = a$Byamiti)

# Fit linear model
Bya_By_CB <- lm(data = a, Byamiti ~ CrocodileBridge)
summary(Bya_By_CB)

# Find the holes
Bya_and_Mal_Wx_holes <- subset(a,is.na(Byamiti))

# Fill the holes
Bya_and_Mal_Wx_holes$Byamiti <- predict(object = Bya_By_CB,newdata = Bya_and_Mal_Wx_holes)

# Check to make sure it looks good
qplot(data = a, x = CrocodileBridge, y = Byamiti) + geom_point(data = Bya_and_Mal_Wx_holes, colour="red", aes(x = CrocodileBridge, y = Byamiti))

rm(a)
rm(Bya_By_CB)

# Houtboschrand Corrections ---------------------------------------------------

# Fill in Houtboschrand. Shimuwini is the closest station.
# Subset first
a <- MakeJustWide(Kruger_Wx_Combined,GoodStation = "Satara", BadStation = "Houtboschrand")

# Test correlation
cor.test(x = a$Satara,y = a$Houtboschrand)

# Fit linear model
model <- lm(data = a, Houtboschrand ~ Satara)
summary(model)

# Find the holes
Hout_and_Sat_Wx_holes <- subset(a,is.na(Houtboschrand))

# Fill the holes
Hout_and_Sat_Wx_holes$Houtboschrand <- predict(object = model,newdata = Hout_and_Sat_Wx_holes)

# Check to make sure it looks good
qplot(data = a, x = Satara, y = Houtboschrand) + geom_point(data = Hout_and_Sat_Wx_holes, colour="red", aes(x = Satara, y = Houtboschrand))

rm(a)
rm(model)

# Pretorioskop Corrections ---------------------------------------------------

# Fill in Pretorioskop. Stolznack is the closest station.
# Subset first
a <- MakeJustWide(Kruger_Wx_Combined,GoodStation = "Stolznek", BadStation = "Pretoriuskop")

# Test correlation
cor.test(x = a$Stolznek,y = a$Pretoriuskop)

# Fit linear model
model <- lm(data = a, Pretoriuskop ~ Stolznek)
summary(model)

# Find the holes
Pret_and_Stol_Wx_holes <- subset(a,is.na(Pretoriuskop))

# Fill the holes
Pret_and_Stol_Wx_holes$Pretoriuskop <- predict(object = model,newdata = Pret_and_Stol_Wx_holes)

# Check to make sure it looks good
qplot(data = a, x = Stolznek, y = Pretoriuskop) + geom_point(data = Pret_and_Stol_Wx_holes, colour="red", aes(x = Stolznek, y = Pretoriuskop))

rm(a)
rm(model)

# Woodlands Corrections ---------------------------------------------------

# Fill in Woodlands. Vlakteplaas is the closest station.
# Subset first
a <- MakeJustWide(Kruger_Wx_Combined,GoodStation = "Vlakteplaas", BadStation = "Woodlands")

# Test correlation
cor.test(x = a$Vlakteplaas,y = a$Woodlands)

# Fit linear model
model <- lm(data = a, Woodlands ~ Vlakteplaas)
summary(model)

# Find the holes
WOO_and_Vla_Wx_holes <- subset(a,is.na(Woodlands))

# Fill the holes
WOO_and_Vla_Wx_holes$Woodlands <- predict(object = model,newdata = WOO_and_Vla_Wx_holes)

# Check to make sure it looks good
qplot(data = a, x = Vlakteplaas, y = Woodlands) + geom_point(data = WOO_and_Vla_Wx_holes, colour="red", aes(x = Vlakteplaas, y = Woodlands))

rm(a)
rm(model)

# Sirheni Corrections ---------------------------------------------------

# Fill in Sirheni. PafuriWenela is the closest station.
# Subset first
a <- MakeJustWide(Kruger_Wx_Combined,GoodStation = "PafuriWenela", BadStation = "Sirheni")

# Test correlation
#cor.test(x = a$PafuriWenela,y = a$Sirheni.)

# Fit linear model
model <- lm(data = a, Sirheni ~ PafuriWenela)
summary(model)

# Find the holes
Sir_and_Paf_Wx_holes <- subset(a,is.na(Sirheni))

# Fill the holes
Sir_and_Paf_Wx_holes$Sirheni <- predict(object = model,newdata = Sir_and_Paf_Wx_holes)

# Check to make sure it looks good
qplot(data = a, x = PafuriWenela, y = Sirheni) + geom_point(data = Sir_and_Paf_Wx_holes, colour="red", aes(x = PafuriWenela, y = Sirheni))

rm(a)
rm(model)


# Combine modeled data ----------------------------------------------------
a <- reshape2::melt(Sir_and_Paf_Wx_holes,id.vars = "Year",variable_name = "Station_Long",value_name = "AnnualPrecip")
b <- reshape2::melt(WOO_and_Vla_Wx_holes,id.vars = "Year",variable_name = "Station_Long",value_name = "AnnualPrecip")
c <- reshape2::melt(Pret_and_Stol_Wx_holes,id.vars = "Year",variable_name = "Station_Long",value_name = "AnnualPrecip")
d <- reshape2::melt(Hout_and_Sat_Wx_holes,id.vars = "Year",variable_name = "Station_Long",value_name = "AnnualPrecip")
e <- reshape2::melt(Bya_and_Mal_Wx_holes,id.vars = "Year",variable_name = "Station_Long",value_name = "AnnualPrecip")
f <- reshape2::melt(Berg_and_Mal_Wx_holes,id.vars = "Year",variable_name = "Station_Long",value_name = "AnnualPrecip")

listDF <- list(a,b,c,d,e,f)

combineModeled <- rbind_all(listDF)
combineModeled$Station_Long <- as.factor(combineModeled$Station_Long)
names(combineModeled)[3] <- "AnnualPrecip"

rm(a,b,c,d,e,f)
rm(Sir_and_Paf_Wx_holes,WOO_and_Vla_Wx_holes,Pret_and_Stol_Wx_holes,Hout_and_Sat_Wx_holes,Bya_and_Mal_Wx_holes,Berg_and_Mal_Wx_holes)

# Merge modeled data into original ----------------------------------------------------

Kruger_Wx_Combined_modeled <- merge(Kruger_Wx_Combined,combineModeled,by=c("Station_Long","AnnualPrecip","Year"),all.x = TRUE)

write.csv(Kruger_Wx_Combined, file = "Weather/Rainfall_LongForm.csv")