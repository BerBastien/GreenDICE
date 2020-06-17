    library("ggExtra")
    library("gridExtra")
    library(wesanderson)
    librar
    y(rpart)
    library("ggridges")
    library("ggplot2")
    library(lubridate)
    library("ggpubr")
    library("tidyverse")
    library(plyr)
    library(randomForest)
    install.packages("readxl")
    library("readxl")


    setwd("C:/Users/bastien/Documents/GitHub/GreenDICE")
    mfp <- read.csv('MFP_naturalcapital.csv')  
    gdp_mfp <- read.csv('GDP_countries_mfp.csv')     
    glimpse(gdp_mfp) 
    glimpse(mfp) 

    years <- factor(gdp_mfp$YearCode)
    years[i]
    for(i in (1:length(levels(years)))){
    gdp_mfp_i <- filter(gdp_mfp, YearCode == years[i])
    trad_mfp_i <- weighted.mean(mfp$traditional_mfp_growth, gdp_mfp_i$AggValue, …)
    mfp_nc_i <- weighted.mean(mfp$mfp_growth_natcap, gdp_mfp_i$AggValue, …)
    tfp_param_i <- mfp_nc_i / trad_mfp_i
    if (i==1){
        tfp_param <-  tfp_param_i
    }else {
        tfp_param <- c(tfp_param, tfp_param_i)
    }
    }
    sd(tfp_param)
    mean(tfp_param)

    #gamma3 = log(tfp_param)/log(K/NC)