################################################################################
## DESCRIPTION: Use meta-regression Bayesian prior tool (MR-BRT) and cascading splines to model prevalence of GNT indicators from 1990 - 2021, and save forecasted SEVs and SDI from 2021 to 2050
## INPUTS: 
#a) Separate prevalence and SEV estimates for each indicator for relevant combinations of location, year, age, sex, and country-level socio-demographic index (SDI)
#b) Separate sets of SEV forecast draws, 500 X each indicator X all relevant combinations of location, year, age, sex, and 
#c) Country-level socio-demographic index (SDI) from 2021 to 2050
#d) Separate population forecasts of children in relevant combinations of location, year, age, sex from 2021 to 2050
## OUTPUTS: 
#a) Cascading spline models for each indicator within the pertinent age groups and sexes, trained on retrospective estimates from 1990 to 2021
#b) SEV forecast data frames (draws) with coresponding country-year SDI forecast present in each row of data ##

## AUTHOR: 
## DATE:
################################################################################


#Getting SEVs for stunting, anemia in WRA, LBW, overweight, EBF, wasting- MOST LIKELY need run in R studio 3.6.1 
#to enable the no random slope on SEV to execute in a cascading spline
library(tidyverse)
library(haven)
library(parallel)
library('reticulate')
#MR-BRT and cascading splines functions require docker. Instructions are available from: https://github.com/ihmeuw-msca/mrtoolr/tree/main
library(reticulate); use_condaenv("mrtool-0.0.1"); library(mrtoolr)
mr <- import("mrtoolr")
library(plyr)
library(dplyr)
library(matrixStats)
library(data.table)
library(zoo)
library(scales)
library(ggplot2)
library(boot)
library(ModelMetrics)


invisible(sapply(list.files("/share/cc_resources/libraries/current/r/", full.names = T), source))
# Directories -------------------------------------------------------------
output_dir <- file.path("FILEPATH")

loc.met <- get_location_metadata(location_set_id = 92, release_id = 9)
#all countries, minus UK, plus the 4 nations of  the UK
loc.met <- loc.met[level==3]
#if we want regional and super regional estimates of SEV as well
#loc.met <- loc.met[level <4]

locs <- unique(loc.met$location_id)
#get the DF in the correct format
prepare_df_format <- function(df, standard.locs){
  
  df[age_group_id == 1, age_group_name := "Under 5"]
  df[age_group_id == 2, age_group_name := "Early Neonatal"]
  df[age_group_id == 3, age_group_name := "Late Neonatal"]
  df[age_group_id == 4, age_group_name := "Post Neonatal"]
  
  #for GBD2019 we need to use agegroupid 5 for overweight in chidren
  df[age_group_id == 5, age_group_name := "1 to 4"]
  df[age_group_id == 388, age_group_name := "1-5 Months"]
  df[age_group_id == 390, age_group_name := "<6 months"]
  df[age_group_id == 389, age_group_name := "6-11 Months"]
  df[age_group_id == 238, age_group_name := "12-23 Months"]
  df[age_group_id == 34, age_group_name := "2-4 Years"]
  df[age_group_id == 24, age_group_name := "15-49 years"]
  df[age_group_id == 164, age_group_name := "Birth"]
  
  df[modelable_entity_id==10556, indicator:="stunting"]
  df[modelable_entity_id==10558, indicator:="wasting"]
  df[modelable_entity_id==20018, indicator:="overweight"]
  df[modelable_entity_id==16282, indicator:="LBW"]
  df[modelable_entity_id==10507, indicator:="anemia"]
  df[modelable_entity_id==20417, indicator:="EBF"]
  
  loc.meta <- standard.locs[, c("location_id", "location_name", "location_type")]
  
  df <- merge(df, loc.meta, by = "location_id")
  
  df[measure_id == 5, measure_name := "Prevalence"]
  df[metric_id == 3, metric_name := "Rate"]
  df[sex_id == 1, sex := "Male"]
  df[sex_id == 2, sex := "Female"]
  df[sex_id == 3, sex := "Both"]
  df[, lyas := paste0(location_id, "_", year_id, "_", age_group_id, "_", sex_id)]
  df[, val := mean(value), by = lyas]
  df[, lower:= quantile(value, probs = .025), by = lyas]
  df[, upper:= quantile(value, probs = .975), by = lyas]
  df <- df[variable == "draw_0"]
  df$variable <- NULL
  df$value <- NULL
  
  return(df)
  
  
  
}
# logit <- function(x){log(x/(1-x))}
# invlogit <- function(x){exp(x)/(1+exp(x))}
id.vars <- c("metric_id", "age_group_id", "location_id", "measure_id", "modelable_entity_id", "sex_id", "year_id", "model_version_id")
hierarchy <- get_location_metadata(92, release_id= 9)

#open retrospective SDI
sdi <- fread("sdi_retrospective.csv")
sdi<-sdi[location_id%in%loc.met$location_id & sex_id==3 & year_id%in%c(1990:2021), c("location_id", "location_name","year_id", "mean_value", "lower_value","upper_value")]
setnames(sdi, old=c("mean","lower", "upper"), new=c("sdi","lower_sdi", "upper_sdi"))

#open retrospective SEVs and prevalence estimates
cgf_df <- readRDS(file.path(output_dir, "cgf_prev2021.RDS"))
sev_df <- readRDS(file.path(output_dir, "SEV_prev2021.RDS"))
overweight_df <- readRDS(file.path(output_dir, "overweight_prev2021.RDS"))


#Forecast models (include all data up to the present:

prepare_forecast_format <- function(df, standard.locs){
  
  df[age_group_id == 1, age_group_name := "Under 5"]
  df[age_group_id == 2, age_group_name := "Early Neonatal"]
  df[age_group_id == 3, age_group_name := "Late Neonatal"]
  df[age_group_id == 4, age_group_name := "Post Neonatal"]
  
  #for GBD2019 we need to use agegroupid 5 for overweight in chidren
  df[age_group_id == 5, age_group_name := "1 to 4"]
  df[age_group_id == 388, age_group_name := "1-5 Months"]
  df[age_group_id == 390, age_group_name := "<6 months"]
  df[age_group_id == 389, age_group_name := "6-11 Months"]
  df[age_group_id == 238, age_group_name := "12-23 Months"]
  df[age_group_id == 34, age_group_name := "2-4 Years"]
  df[age_group_id == 24, age_group_name := "15-49 years"]
  df[age_group_id == 164, age_group_name := "Birth"]
  
  loc.meta <- standard.locs[, c("location_id", "location_name", "location_type")]
  
  df <- merge(df, loc.meta, by = "location_id")
  
  df[sex_id == 1, sex := "Male"]
  df[sex_id == 2, sex := "Female"]
  df[sex_id == 3, sex := "Both"]
  df[, lyas := paste0(location_id, "_", year_id, "_", age_group_id, "_", sex_id, "_", scenario)]
  df[, sev := mean(value), by = lyas]
  df[, lower_sev:= quantile(value, probs = .025), by = lyas]
  df[, upper_sev:= quantile(value, probs = .975), by = lyas]
  df <- df[draw == 0]
  df$variable <- NULL
  df$draw <- NULL
  
  return(df)
}

prepare_forecast_format_draws <- function(df, standard.locs){
  
  df[age_group_id == 1, age_group_name := "Under 5"]
  df[age_group_id == 2, age_group_name := "Early Neonatal"]
  df[age_group_id == 3, age_group_name := "Late Neonatal"]
  df[age_group_id == 4, age_group_name := "Post Neonatal"]
  
  #for GBD2019 we need to use agegroupid 5 for overweight in chidren
  #fore GBD2021 we need to use agegroupid 34 for overweight in chidren
  df[age_group_id == 5, age_group_name := "1 to 4"]
  df[age_group_id == 388, age_group_name := "1-5 Months"]
  df[age_group_id == 390, age_group_name := "<6 months"]
  df[age_group_id == 389, age_group_name := "6-11 Months"]
  df[age_group_id == 238, age_group_name := "12-23 Months"]
  df[age_group_id == 34, age_group_name := "2-4 Years"]
  df[age_group_id == 24, age_group_name := "15-49 years"]
  df[age_group_id == 164, age_group_name := "Birth"]
  
  loc.meta <- standard.locs[, c("location_id", "location_name", "location_type")]
  
  df <- merge(df, loc.meta, by = "location_id")
  
  df[sex_id == 1, sex := "Male"]
  df[sex_id == 2, sex := "Female"]
  df[sex_id == 3, sex := "Both"]
  # df[, lyas := paste0(location_id, "_", year_id, "_", age_group_id, "_", sex_id, "_", scenario)]
  # df[, sev := mean(value), by = lyas]
  # df[, lower_sev:= quantile(value, probs = .025), by = lyas]
  # df[, upper_sev:= quantile(value, probs = .975), by = lyas]
  # df <- df[draw == 0]
  # df$variable <- NULL
  # df$draw <- NULL
  
  return(df)
}

pop_forecast<- fread("population_forecast.csv")
pop_forecast <- pop_forecast[location_id %in% locs,]
setnames(pop_forecast, old="value", new= "population")

pop_forecast<- pop_forecast[order(location_id, age_group_id, sex_id, year_id)]


# SDI forecasts
sdi_forecast <- fread("sdi_forecast.csv")
sdi_forecast<- sdi_forecast[location_id%in%locs & year_id%in%c(2021:2050) & scenario== 0,]
setnames(sdi_forecast, "mean", "sdi") 

#stunting forecasts:
stunt_forecast_draws<- fread("FILEPATH")
stunt_forecast_draws<- stunt_forecast_draws[year_id%in%c(2021:2050) & draw< 500 & scenario == 0,]
stunt_sev_forecast_draws <- prepare_forecast_format_draws(df = stunt_forecast_draws, loc.met[level==3])
rm(stunt_forecast_draws)

stunt_sev_forecast_draws <- merge(stunt_sev_forecast_draws, hierarchy[, .(location_id, super_region_name, region_name)], all.x=T, by = "location_id")
setnames(stunt_sev_forecast_draws, old= "value" , new= "sev")
stunt_sev_forecast_draws <- stunt_sev_forecast_draws[location_id%in%loc.met$location_id, c("location_id", "year_id", "age_group_id","sex_id","scenario","super_region_name","region_name","sev", "draw")]
#saveRDS(stunt_sev_forecast_draws, file.path(paste0("FILEPATH")))
stunt_sev_forecast_draws<- readRDS("FILEPATH")

#mean and UIs:
stunt_sev_forecast <- prepare_forecast_format(df = stunt_forecast_draws, loc.met[level==3])

#columns to keep from retrospective dataframe:
keep_cols<- c("location_id", "location_name", "sex_id", "rei_id", "age_group_id", "super_region_name", "region_name", "age_group_name")
#create stunt_forecast_ref for agegroups 2 and 3, using forecasted sevs from agegroup4
stunt_sev_forecast_23 <- merge(cgf_df[rei_id== 241 & year_id==1990 & age_group_id%in%c(2,3), ..keep_cols], stunt_sev_forecast_draws[age_group_id==4, -c("age_group_id")], by=c("location_id", "sex_id", "region_name", "super_region_name"), allow.cartesian = TRUE)
#saveRDS(stunt_sev_forecast_23, file.path(paste0("FILEPATH")))
stunt_sev_forecast_23 <- readRDS("/ihme/mnch/cgf/model/iterative/papers/gnt_paper/Forecasting/Forecasts/sev/stunt_sev_forecast_23_20240119.RDS")

#merge key columns from stunt_df, but just from one of the years.
stunt_sev_forecast_draws <- merge(cgf_df[rei_id== 241 & year_id==1990, ..keep_cols], stunt_sev_forecast_draws, by=c("location_id", "age_group_id", "sex_id", "region_name", "super_region_name"))

#rbind the artificial sev forecasts for agegroups 2 and 3 with the real sev forecasts for groups 4, 5:
stunt_sev_forecast <- rbind(stunt_sev_forecast_23, stunt_sev_forecast_draws)
rm(stunt_sev_forecast_23, stunt_sev_forecast_draws)
#merge sdi and sev forecast DFs- use reference scenario population
stunt_sev_forecast <- merge(stunt_sev_forecast, sdi_forecast[, c("location_id", "year_id", "sdi", "scenario")], by=c("location_id", "year_id", "scenario"))
#merge population forecasts
stunt_sev_forecast <- merge(stunt_sev_forecast, pop_forecast, by=c("location_id", "year_id", "age_group_id", "sex_id"))


#create stunting forecasts for all 4 sex/agegroup combinations
sid <- 1
agid <- 2
scenid <- 1

saveRDS(stunt_sev_forecast, file.path(paste0("FILEPATH")))

#wasting forecasts
waste_forecast_draws<- fread("FILEPATH")
waste_forecast_draws<- waste_forecast_draws[year_id%in%c(2021:2050) & draw< 500 & scenario == 0,]

waste_sev_forecast_draws <- prepare_forecast_format_draws(df = waste_forecast_draws, loc.met[level==3])
waste_sev_forecast_draws <- merge(waste_sev_forecast_draws, hierarchy[, .(location_id, super_region_name, region_name)], all.x=T, by = "location_id")
setnames(waste_sev_forecast_draws, old= "value" , new= "sev")
waste_sev_forecast_draws <- waste_sev_forecast_draws[location_id%in%loc.met$location_id, c("location_id", "year_id", "age_group_id","sex_id","scenario","super_region_name","region_name","sev", "draw")]

#create waste_forecast_ref for agegroups 2 and 3, using forecasted sevs from agegroup4
waste_sev_forecast_23 <- merge(cgf_df[year_id==1990 & age_group_id%in%c(2,3) & rei_id==240, ..keep_cols], waste_sev_forecast_draws[age_group_id==4, -c("age_group_id")], by=c("location_id", "sex_id", "region_name", "super_region_name"), allow.cartesian = TRUE)


#merge key columns from cgf_df, but just from one of the years.
waste_sev_forecast_draws <- merge(cgf_df[year_id==1990  & rei_id==240, ..keep_cols], waste_sev_forecast_draws, by=c("location_id", "age_group_id", "sex_id", "region_name", "super_region_name"))

#rbind the artificial sev forecasts for agegroups 2 and 3 with the reali sev forecasts for groups 4, 5:
waste_sev_forecast <- rbind(waste_sev_forecast_23, waste_sev_forecast_draws)

#merge sdi and sev forecast DFs- match on scenario
waste_sev_forecast <- merge(waste_sev_forecast, sdi_forecast, by=c("location_id", "year_id", "scenario"))
#merge population forecasts
waste_sev_forecast <- merge(waste_sev_forecast, pop_forecast, by=c("location_id", "year_id", "age_group_id", "sex_id"))
rm(waste_sev_forecast_23, waste_sev_forecast_draws, waste_forecast_draws)

saveRDS(waste_sev_forecast, file.path(paste0("FILEPATH")))
waste_sev_forecast <- readRDS(paste0("FILEPATH"))


#LBW forecasts:
LBW_forecast_draws<- fread("FILEPATH")
LBW_forecast_draws<- LBW_forecast_draws[year_id%in%c(2021:2050) & draw< 500 & scenario == 0,]

LBW_sev_forecast_draws <- prepare_forecast_format_draws(df = LBW_forecast_draws, loc.met[level==3])
LBW_sev_forecast_draws <- merge(LBW_sev_forecast_draws, hierarchy[, .(location_id, super_region_name, region_name)], all.x=T, by = "location_id")
setnames(LBW_sev_forecast_draws, old= "value" , new= "sev")
LBW_sev_forecast_draws <- LBW_sev_forecast_draws[location_id%in%loc.met$location_id, c("location_id", "year_id", "age_group_id","sex_id","scenario","super_region_name","region_name","sev", "draw")]

#merge key columns from sev_df, but just from one of the years.
sev_df <- fread("GNT_SEVs.csv")
setnames(sev_df, c("val", "upper", "lower"), c("sev", "upper_sev", "lower_sev"))
prev_df <- fread("GNT_Prev.csv")
setnames(prev_df, c("val", "upper", "lower"), c("prevalence", "upper_prev", "lower_prev"))
sev_df= merge(sev_df[, ], prev_df[, c("location_id", "year_id", "age_group_id","sex_id","prevalence", "upper_prev", "lower_prev")], by= c("location_id", "year_id", "age_group_id", "sex_id"))

LBW_sev_forecast_draws <- merge(sev_df[year_id==1990  & rei_id==335,
                                       c("location_id","location_name", "sex_id", "measure_id","metric_id", "rei_id", 
                                         "age_group_name", "age_group_id", "super_region_name", "region_name")], 
                                LBW_sev_forecast_draws[age_group_id==2,], by=c("location_id", "age_group_id", "sex_id", "region_name", "super_region_name"))

#merge sdi and sev forecast DFs- match on scenario
LBW_sev_forecast_draws <- merge(LBW_sev_forecast_draws, sdi_forecast, by=c("location_id", "year_id", "scenario"))

#change the age_id for LBW prediction, to specify group for prevalence estimates
LBW_sev_forecast_draws[, age_group_id:= 164]
LBW_sev_forecast_draws[, age_group_name:= "Birth"]

#load forecasts and retrospective estimates of livebirths
livebirth_forecast<- fread("births_forecast.csv" )
births_retrospective <- fread("births_retrospective.csv")
livebirth_forecast<- rbind(births_retrospective[sex_id!=3,-"run_id"], livebirth_forecast, fill= TRUE)

LBW_sev_forecast <- merge(LBW_sev_forecast_draws, livebirth_forecast[, -c("age_group_name", "scenario")], by=c("location_id", "year_id", "sex_id", "age_group_id"))
rm(LBW_sev_forecast_draws, LBW_forecast_draws)

saveRDS(LBW_sev_forecast, file.path(paste0("FILEPATH")))
LBW_sev_forecast <- readRDS(paste0("FILEPATH"))


#EBF forecasts:
EBF_forecast_draws<- fread("FILEPATH")
EBF_forecast_draws<- EBF_forecast_draws[year_id%in%c(2021:2050) & draw< 500 & scenario == 0,]

EBF_sev_forecast_draws <- prepare_forecast_format_draws(df = EBF_forecast_draws, loc.met[level==3])
EBF_sev_forecast_draws <- merge(EBF_sev_forecast_draws, hierarchy[, .(location_id, super_region_name, region_name)], all.x=T, by = "location_id")
setnames(EBF_sev_forecast_draws, old= "value" , new= "sev")
EBF_sev_forecast_draws <- EBF_sev_forecast_draws[location_id%in%loc.met$location_id, c("location_id", "year_id", "age_group_id","sex_id","scenario","super_region_name","region_name","sev", "draw")]
EBF_sev_forecast_draws <- EBF_sev_forecast_draws[age_group_id==3,]

#merge key columns from sev_df, but just from one of the years.
EBF_sev_forecast <- merge(sev_df[year_id==1990  & rei_id==136,
                                 c("location_id","location_name", "sex_id", "measure_id","metric_id", "rei_id", 
                                   "age_group_name", "age_group_id", "super_region_name", "region_name")], 
                          EBF_sev_forecast_draws, by=c("location_id", "age_group_id", "sex_id", "region_name", "super_region_name"))


#merge sdi and sev forecast DFs- match on scenario
EBF_sev_forecast <- merge(EBF_sev_forecast, sdi_forecast, by=c("location_id", "year_id", "scenario"))
#merge population forecasts
EBF_sev_forecast <- merge(EBF_sev_forecast, pop_forecast, by=c("location_id", "year_id", "age_group_id", "sex_id"))
rm(EBF_sev_forecast_draws, EBF_forecast_draws)

#saveRDS(EBF_sev_forecast, file.path(paste0("FILEPATH")))
EBF_sev_forecast <- readRDS(paste0("FILEPATH"))

#overweight forecasts:
adult_bmi_forecast_draws <- fread("FILEPATH")
adult_bmi_forecast_draws<- adult_bmi_forecast_draws[year_id%in%c(2021:2050) & draw< 500 & scenario == 0,]
#merge population so we can take the population weighted mean of the SEV from age groups 20-49
adult_bmi_forecast_draws<- merge(adult_bmi_forecast_draws[age_group_id%in%c(9:14) & scenario==0,], pop_forecast, by= c("location_id", "sex_id", "year_id", "age_group_id"))

adult_bmi_forecast_draws[, sev_adults:= weighted.mean(value, population), by= c("location_id", "sex_id", "year_id", "draw")]
adult_bmi_forecast_draws<- adult_bmi_forecast_draws[age_group_id==9,]
adult_bmi_forecast_draws[, age_group_id:= 34]
adult_bmi_forecast_draws[, age_group_name := "2-4 Years"]

adult_bmi_forecast_DF <- prepare_forecast_format_draws(df = adult_bmi_forecast_draws, loc.met[level==3])
adult_bmi_forecast_DF <- merge(adult_bmi_forecast_DF, hierarchy[, .(location_id, super_region_name, region_name)], all.x=T, by = "location_id")
adult_bmi_forecast_DF <- adult_bmi_forecast_DF[location_id%in%loc.met$location_id, c("location_id", "year_id", "age_group_id","sex_id","scenario","super_region_name","region_name","sev_adults", "draw")]

#add child SEVs since these have now been forecasted...
child_overweight_forecast_draws<- fread("FILEPATH")
child_overweight_forecast_draws<- child_overweight_forecast_draws[year_id%in%c(2021:2050) & draw< 500 & scenario == 0,]
child_overweight_sev_forecast_draws <- prepare_forecast_format_draws(df = child_overweight_forecast_draws, loc.met[level==3])
setnames(child_overweight_sev_forecast_draws, old= "value" , new= "sev_kids")

#merge child overweight SEVs with the adult overweight SEVs that have already been formatted:
overweight_forecast_DF<- merge(adult_bmi_forecast_DF, child_overweight_sev_forecast_draws[age_group_id==6,c("location_id", "draw", "scenario", "sex_id", "year_id", "sev_kids")], by= c("location_id", "sex_id", "year_id", "draw", "scenario") )

#merge sdi and sev forecast DFs- match on scenario
overweight_forecast_DF <- merge(overweight_forecast_DF, sdi_forecast, by=c("location_id", "year_id", "scenario"))

overweight_forecast_DF[,rei_id:=371]
#saveRDS(overweight_forecast_DF, file.path(paste0("FILEPATH")))
overweight_forecast_DF<- readRDS(paste0("FILEPATH"))



#anemia (iron deficiency SEVs) forecasts:
anem_forecast_draws<- fread("FILEPATH")
anem_forecast_draws<- anem_forecast_draws[year_id%in%c(2021:2050) & draw< 500 & scenario == 0,]

anem_sev_forecast_draws <- prepare_forecast_format_draws(df = anem_forecast_draws, loc.met[level==3])
anem_sev_forecast_draws <- merge(anem_sev_forecast_draws, hierarchy[, .(location_id, super_region_name, region_name)], all.x=T, by = "location_id")
setnames(anem_sev_forecast_draws, old= "value" , new= "sev")
anem_sev_forecast <- anem_sev_forecast_draws[location_id%in%loc.met$location_id, c("location_id", "year_id", "age_group_id","sex_id","scenario","super_region_name","region_name","sev", "draw")]

#merge sdi and sev forecast DFs- match on scenario
anem_sev_forecast <- merge(anem_sev_forecast, sdi_forecast, by=c("location_id", "year_id", "scenario"))
#merge population forecasts
anem_sev_forecast <- merge(anem_sev_forecast, pop_forecast, by=c("location_id", "year_id", "age_group_id", "sex_id"))
rm(anem_sev_forecast_draws, anem_forecast_draws)

anem_sev_forecast<- anem_sev_forecast[age_group_id%in%c(8:14),]
anem_sev_forecast[,rei_id:=95]
saveRDS(anem_sev_forecast, file.path(paste0("FILEPATH")))






#Stunting and Wasting (CGF) full cascade spline runs to be saved before running array scripts (sex and age-specific)
#Run the cascading spline models on 1990-2021 data so they may be opened by the array jobs for prediction/forecasting
library(parallel)
cl <- makeCluster(getOption("cl.cores", 4))
lapply(c(240,241), function(r){
  lapply(c(2:5), function(agid){  
    library(data.table)
    #MR-BRT and cascading splines functions require docker. Instructions are available from: https://github.com/ihmeuw-msca/mrtoolr/tree/main
    library(reticulate); use_condaenv("mrtool-0.0.1"); library(mrtoolr)
    mr <- import("mrtoolr")
    library(mrbrt001)
    library(boot)
    library(plyr)
    library(dplyr)


    standard <- "cgf_full"
    cascade <- "cgf_cascade_full"
    
    cgf_df <- readRDS(file.path("FILEPATH/cgf_prev2021.RDS"))
    output_dir2<-"FILEPATH"
    fit_global_full <- function(model_data, r, agid, sid) {
      #subset to only contain both sexes and split into training and testing dataframes
      train_df<- model_data[rei_id==r & sex_id== sid & age_group_id== agid,]

      #looking at sev distribution
      hist(model_data[rei_id==r & sex_id== sid & age_group_id== agid,]$sev)
      #training data
      dat_loc <- MRData()
      dat_loc$load_df(
        data = train_df,
        col_obs = "logit_prev", col_obs_se = "se",
        col_covs = list("sdi", "sev"), col_study_id = "location_id"
      )
      
      #"standard MRBRT model
      stand_mod <- MRBRT(
        data = dat_loc,
        cov_models = list(
          LinearCovModel("intercept", use_re = TRUE),
          LinearCovModel("sdi", use_re = FALSE),
          LinearCovModel("sev", use_re = FALSE, 
                         use_spline = TRUE,
                         #fit spline to every quartile of the SEV distribution
                         spline_knots = array(seq(0, 1, by = 0.25)), # this tells is to put a spline every 0.25
                         spline_degree = 3L,
                         spline_knots_type = 'frequency',# this specifies to put the knots literally along the domain. Other option is frequency for it to be data density driven
                         spline_r_linear = TRUE,
                         spline_l_linear = FALSE
          )
        ))
      stand_mod$fit_model(inner_print_level = 5L, inner_max_iter = 1000L,  outer_max_iter= 500L)
      estimated_beta <- stand_mod$summary()[[1]][1, "sdi"]
      stand_mod$cov_models[[which(stand_mod$cov_names == "sdi")]]$prior_beta_uniform <- matrix(rep(estimated_beta, 2), ncol = 1)
      
      return(stand_mod)
      
    }
    fit_cascade_full <- function(mod_global, model_data, output_dir, r, agid, sid) {
      
      train_df<- model_data[rei_id==r & sex_id== sid & age_group_id== agid,]
      
      model_label_tmp <- paste0("cascade_prevalence", unique(train_df[rei_id==r,]$rei_id), "_", unique(train_df[rei_id==r,]$age_group_id), "_", unique(train_df$sex))
      
      thetas <- c(3,3)
      cascade_fit <- run_spline_cascade(
        stage1_model_object = mod_global, 
        df = train_df, 
        col_obs = "logit_prev",
        col_obs_se = "se", 
        col_study_id = "location_id", 
        stage_id_vars = c("region_name", "location_id"),
        thetas = thetas,
        gaussian_prior = TRUE,
        output_dir = output_dir,
        model_label = model_label_tmp, 
        overwrite_previous = TRUE
      )
      return(cascade_fit)
    }
    
    for(sid in c(1,2)){
      
      #fit models
      cgf_full <- fit_global_full(cgf_df,r, agid, sid)
      cgf_cascade_full <- fit_cascade_full(cgf_full,cgf_df,output_dir2, r, agid,sid)
      
    }
  })
})

#EBF, LBW, and anaemia
#Run the cascading spline models on 1990-2021 data so they may be opened by the array jobs for prediction/forecasting
for (r in c(136, 335, 95)){
  
  if (r==335) {
    ind<-"lbw"
    age_list <- 164
    sex_list<- c(1,2)
    standard <- "mod_full"
    cascade <- "mod_cascade_full"
    thetas <- c(3,3)
    sev_df <- readRDS(file.path("FILEPATH/SEV_prev2021.RDS"))
    #change the age_group_id to match LBW's 164
    sev_df[rei_id==335 ,age_group_id:=164]
  }
  if (r==136) {
    ind<- "ebf"
    age_list<- 3
    sex_list<- c(1,2)
    standard <- "mod_full"
    cascade <- "mod_cascade_full"
    thetas <- c(3,3)
    sev_df <- readRDS(file.path("FILEPATH/SEV_prev2021.RDS"))
  }
  
  if (r==95) {
    ind<- "anem"
    age_list<- c(8:14)
    sex_list<- 2
    standard <- "mod_full"
    cascade <- "mod_cascade_full"
    thetas <- c(5,3)
    sev_df <- readRDS(file.path("FILEPATH/anem_prev2021.RDS"))
    
  }
  
  for(agid in age_list){
    
    for(sid in sex_list){ 
      
      output_dir2<-"FILEPATH"
      fit_global_full <- function(model_data, r, agid, sid) {
        #subset to only contain both sexes and split into training and testing dataframes
        train_df<- model_data[rei_id==r & sex_id== sid & age_group_id== agid,]
        #looking at sev distribution
        hist(model_data[rei_id==r & sex_id== sid & age_group_id== agid,]$sev)
        #training data
        dat_loc <- MRData()
        dat_loc$load_df(
          data = train_df,
          col_obs = "logit_prev", col_obs_se = "se",
          col_covs = list("sdi", "sev"), col_study_id = "location_id"
        )
     
        #"standard MRBRT model
        stand_mod <- MRBRT(
          data = dat_loc,
          cov_models = list(
            LinearCovModel("intercept", use_re = TRUE),
            LinearCovModel("sdi", use_re = FALSE),
            LinearCovModel("sev", use_re = FALSE, 
                           use_spline = TRUE,
                           #fit spline to every quartile of the SEV distribution
                           spline_knots = array(seq(0, 1, by = 0.25)), # this tells is to put a spline every 0.25
                           spline_degree = 3L,
                           spline_knots_type = 'frequency',# this specifies to put the knots literally along the domain. Other option is frequency for it to be data density driven
                           spline_r_linear = TRUE,
                           spline_l_linear = FALSE
            )
          ))
        stand_mod$fit_model(inner_print_level = 2L, inner_max_iter = 2000L,  outer_max_iter= 1000L)
        estimated_beta <- stand_mod$summary()[[1]][1, "sdi"]
        stand_mod$cov_models[[which(stand_mod$cov_names == "sdi")]]$prior_beta_uniform <- matrix(rep(estimated_beta, 2), ncol = 1)
        
        return(stand_mod)
        
      }
      
      
      fit_cascade_full <- function(mod_global, model_data, output_dir, r, agid, sid) {
        
        train_df<- model_data[rei_id==r & sex_id== sid & age_group_id== agid,]
        
        
        model_label_tmp <- paste0("cascade_prevalence", unique(train_df[rei_id==r,]$rei_id), "_", unique(train_df[rei_id==r,]$age_group_id), "_", unique(train_df$sex))
        
        thetas <- thetas
        cascade_fit <- run_spline_cascade(
          stage1_model_object = mod_global, 
          df = train_df, 
          col_obs = "logit_prev",
          col_obs_se = "se", 
          col_study_id = "location_id", 
          stage_id_vars = c("region_name", "location_id"),
          thetas = thetas,
          gaussian_prior = TRUE,
          output_dir = output_dir,
          model_label = model_label_tmp, 
          overwrite_previous = TRUE
        )
        return(cascade_fit)
      }
      
      
      mod_full <- fit_global_full(sev_df,r, agid, sid)
      sex <- ifelse(sid==1, "Male", "Female")
      mod_cascade_full <- fit_cascade_full(mod_full,sev_df,output_dir2, r, agid,sid)
    }
  }
}
#output from array jobs will need to be intercept shifted before aggregation by sex.



#Overweight
#Run the cascading spline models on 1990-2021 data so they may be opened by the array jobs for prediction/forecasting
for(sid in c(1,2)){
  r<- 371
  agid <- 34
  ind<-"overweight"
  standard <- "overweight_full"
  cascade <- "overweight_cascade_full"
  
  overweight_df <- readRDS(file.path("FILEPATH"))#overweight prevalence for the year 2021
  
  output_dir2<-"FILEPATH"
  #fit standard MR-BRT model
  fit_global_full <- function(model_data, r, agid, sid) {
    train_df<- model_data[rei_id==r & sex_id== sid & age_group_id== agid,]
    #looking at mean adult BMI distribution
    hist(model_data[rei_id==r & sex_id== sid & age_group_id== agid,]$sev_kids)
    #training data
    dat_loc <- MRData()
    dat_loc$load_df(
      data = train_df,
      col_obs = "logit_prev", col_obs_se = "se",
      col_covs = list("sdi", "sev_kids"), col_study_id = "location_id"
    )
    
    #"standard MRBRT model
    stand_mod <- MRBRT(
      data = dat_loc,
      cov_models = list(
        LinearCovModel("intercept", use_re = TRUE),
        LinearCovModel("sdi", use_re = FALSE),
        LinearCovModel("sev_kids", use_re = FALSE, 
                       use_spline = TRUE,
                       #fit spline to every quartile of the mean BMI distribution
                       spline_knots = array(seq(0, 1, by = 0.25)), # this tells is to put a spline every 0.25
                       spline_degree = 3L,
                       spline_knots_type = 'frequency',# this specifies to put the knots literally along the domain. Other option is frequency for it to be data density driven
                       spline_r_linear = TRUE,
                       spline_l_linear = FALSE
        )
      ))
    stand_mod$fit_model(inner_print_level = 5L, inner_max_iter = 1000L,  outer_max_iter= 500L)
    
    return(stand_mod)
    
  }
  #fit cascading spline model
  fit_cascade_full <- function(mod_global, model_data, output_dir, r, agid, sid) {
    
    train_df<- model_data[rei_id==r & sex_id== sid & age_group_id== agid,]
    
    
    model_label_tmp <- paste0("cascade_prevalence", unique(train_df[rei_id==r,]$rei_id), "_", unique(train_df[rei_id==r,]$age_group_id), "_", unique(train_df$sex))
    thetas <- c(1,3)
    cascade_fit <- run_spline_cascade(
      stage1_model_object = mod_global, 
      df = train_df, 
      col_obs = "logit_prev",
      col_obs_se = "se", 
      col_study_id = "location_id", 
      stage_id_vars = c("region_name", "location_id"),
      thetas = thetas,
      gaussian_prior = TRUE,
      output_dir = output_dir,
      model_label = model_label_tmp, 
      overwrite_previous = TRUE
    )
    return(cascade_fit)
  }   
  
  #run model on sex-specific overweight and agegroup
  mod_full <- fit_global_full(overweight_df,r, agid, sid)
  sex <- ifelse(sid==1, "Male", "Female")
  mod_cascade_full <- fit_cascade_full(mod_full,overweight_df,output_dir2, r, agid,sid)
  
}