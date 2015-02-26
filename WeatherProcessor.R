library(reshape)
library(plyr)
library(ggplot2)
library(lattice)

setwd("/Users/danielgodwin/Dropbox/Graduate School/Dissertation/Chapter 3 - Growth/TreeRings/")

sitesWx <- read.csv("Weather/Kruger_Stations1980_2008.csv")

sitesWx_long <- melt(sitesWx,id.vars="Year",variable_name="Station")
names(sitesWx_long)[3] <- "AnnualPrecip"

newWx <- read.csv("Weather/monthly_rainfall_2006_2010_paf-moo-sat-pre.txt")
newWx_annual <- ddply(newWx,.(YEAR,STATION),summarise,AnnualPrecip = sum(SumOfMM))
names(newWx_annual) <- c("Year","Station","AnnualPrecip")

newWx_annual_subset <- subset(newWx_annual,Year >= 2009)

Kruger_Wx_Combined <- rbind(sitesWx_long,newWx_annual_subset)

lookupStation <- data.frame(Station = unique(Kruger_Wx_Combined$Station),Station_Long = c("Satara","Berg En Dal","Byamiti","Houtboschrand","Kingfisherspruit","Crocodile Bridge","Letaba","Mahlangeni","Malelane","Mooiplaas","Mopani","Nwanetsi","Olifants",
                                                                                          "Lower Sabie","Pafuri","Pafuri (Wenela)","Phalaborwa","Pretoriuskop","Punda Maria","Shangoni","Shingwedzi","Shimuwini","Sirheni","Skukuza","Stolznek","Talamati",
                                                                                          "Tshokwane","Vlakteplaas","Woodlands"))

Kruger_Wx_Combined <- merge(Kruger_Wx_Combined,lookupStation,by="Station")


WRF_Wx <- read.csv("Weather/WRF_rainfall_annual_totals.csv")
names(WRF_Wx)[2] <- "AnnualPrecip"

WRF_Wx$Station <- "WRF"
WRF_Wx$Station_Long <- "Wits Rural Facility"

Kruger_Wx_Combined <- rbind(Kruger_Wx_Combined,WRF_Wx)
