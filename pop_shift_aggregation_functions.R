################################################################################
## DESCRIPTION: Intercept-shift the forecasted population estimates (age groups 34 and 388) from 2021 onwards by taking the natural-space differences 
#between the year 2021 forecast and the 2021 GBD estimate. Then aggregate the shifted population forecasts and livebirth forecasts to the desired location, 
#sex, age groups. Merge with the previously shifted population forecasts and sum component populations to create the new age group 390 (under 6 months).

## INPUTS: 
#a) Data frame containing country-age-sex-year specific mean population forecasts from 2022 to 2050 for most age groups relevant to GNT indicators (except for the new GBD 2021 age groups).
#b) Data frame containing sex-specific livebirth forecasts
#c) Separate data frames containing forecasted country-age-sex-year mean population from 2021 to 2050 for age groups 34 (2-4 years) and 388 (1-5 months) [required to create new relevant age groups]
#d) GBD 2021 country-age-sex specific mean population estimates for year 2021 in age group ids 164, 1, 2,3,4, 5, 388, 8:14, 24, 34

## OUTPUTS: 
#a) Data frame containing country-age-sex-year specific mean population forecasts from 2021 to 2050 for all age groups relevant to GNT indicators 

invisible(sapply(list.files("FILEPATH", full.names = T), source))

library(matrixStats)
library(data.table)
library(zoo)
library(plyr)
library(dplyr)
require(boot)

loc.met <- get_location_metadata(location_set_id = 91, release_id = 9)
loc.met <- loc.met[level %in% c(0,1, 2, 3)]
locs <- loc.met$location_id

id.vars <- c("metric_id", "age_group_id", "location_id", "measure_id", "modelable_entity_id", "sex_id", "year_id", "model_version_id")
hierarchy <- get_location_metadata(92, release_id= 9) #location hierarchy
pop_2021 <- get_population(location_id = locs, release_id = 9, location_set_id = 91, year_id = 2021, age_group_id = c(164, 1, 2,3,4, 5, 388, 8:14, 24, 34), sex_id = c(1,2,3))#get populations
setnames(pop_2021, "population", "value")
#create GBD 2021 population estimate for age group 390 by summing across component age groups
pop_390_2021<- pop_2021[age_group_id%in%c(388, 2, 3)]
pop_390_2021[, pop_tot:= sum(value), by= c("location_id", "sex_id", "year_id")]
pop_390_2021<- pop_390_2021[age_group_id==388, -"value"]
pop_390_2021[, age_group_id:= 390]
setnames(pop_390_2021, "pop_tot", "value")
pop_2021<- rbind(pop_2021, pop_390_2021)



pop_forecast<- fread("FILEPATH/population.csv")
pop_forecast<- pop_forecast[year_id>2020 & scenario==0 & age_group_id %in% c(164, 1, 2,3,4, 5, 8:14, 24) & location_id%in%locs,]
setnames(pop_forecast, "mean", "value")
summary(pop_forecast$value)
#merge y2021 population with forecasts:
pop_forecast<- rbind(pop_2021[,-"run_id"], pop_forecast[age_group_id %in% c(164, 1, 2,3,4, 5, 8:14, 24), ], fill= TRUE)

#function to provide the age group name, given the age_group_id
ag_name_format<- function(age_group_id){
  if(age_group_id == 1){age_group_name<- "Under 5"}
  if(age_group_id == 2){age_group_name<- "Early Neonatal"}
  if(age_group_id == 3){age_group_name<- "Late Neonatal"}
  if(age_group_id == 4){age_group_name<- "Post Neonatal"}
  if(age_group_id == 5){age_group_name<- "1 to 4"}
  if(age_group_id == 6){age_group_name<- "5 to 9"}
  if(age_group_id == 7){age_group_name<- "10 to 14"}
  if(age_group_id == 8){age_group_name<- "15 to 19"}
  if(age_group_id == 9){age_group_name<- "20 to 24"}
  if(age_group_id == 10){age_group_name<- "25 to 29"}
  if(age_group_id == 11){age_group_name<- "30 to 34"}
  if(age_group_id == 12){age_group_name<- "35 to 39"}
  if(age_group_id == 13){age_group_name<- "40 to 44"}
  if(age_group_id == 14){age_group_name<- "45 to 49"}
  if(age_group_id == 15){age_group_name<- "50 to 54"}
  if(age_group_id == 16){age_group_name<- "55 to 59"}
  if(age_group_id == 17){age_group_name<- "60 to 64"}
  if(age_group_id == 18){age_group_name<- "65 to 69"}
  if(age_group_id == 19){age_group_name<- "70 to 74"}
  if(age_group_id == 20){age_group_name<- "75 to 79"}
  if(age_group_id == 30){age_group_name<- "80 to 84"}
  if(age_group_id == 31){age_group_name<- "85 to 89"}
  if(age_group_id == 32){age_group_name<- "90 to 94"}
  if(age_group_id == 235){age_group_name<- "95 plus"}
  if(age_group_id == 257){age_group_name<- "20-79 years"}
  if(age_group_id == 388){age_group_name<- "1-5 Months"}
  if(age_group_id == 390){age_group_name<- "<6 months"}
  if(age_group_id == 389){age_group_name<- "6-11 Months"}
  if(age_group_id == 238){age_group_name<- "12-23 Months"}
  if(age_group_id == 34){age_group_name<- "2-4 Years"}
  if(age_group_id == 24){age_group_name<- "15-49 years"}
  if(age_group_id == 164){age_group_name<- "Birth"}
  return(age_group_name)
}

#function to intercept shift the population forecast
pop_intercept_shift<- function(pop_forecast){
  forecast<- pop_forecast[year_id>2020 & scenario==0 & age_group_id %in% c(164, 1, 2,3,4, 5, 388, 8:14, 24, 34),]
  gbd_estimate_2021 <- get_population(location_id = locs, release_id = 9, location_set_id = 91, year_id = 2021, age_group_id = c(164, 1, 2,3,4, 5, 388, 8:14, 24, 34), sex_id = c(1,2,3), run_id = 359)
  #determine the difference between the mean population forecast for Y2021 and the mean population estimate from GBD 2021
  temp<- merge(forecast, gbd_estimate_2021[, c("location_id","year_id","age_group_id", "sex_id","population")], by=c("location_id","year_id","age_group_id", "sex_id"))
  #additive intercept shift
  temp[year_id==2021, shift_difference:= population- value] 
  #multiplicative intercept shift
  temp[year_id==2021, shift_relative:= population/value] 
  
  
  forecast<- merge(forecast, temp[,c("location_id","age_group_id", "sex_id","shift_difference", "shift_relative")], by= c("location_id","age_group_id", "sex_id"))
  
  forecast[, population:= value + shift_difference]
  #forecast[, population:= value* shift_relative]
  forecast<- forecast[,-c("value", "shift_difference", "shift_relative")]
  setnames(forecast, "population", "value")
  return(forecast)
  
}

#Previously shifted birth forecasts:
birth_forecast<- fread("FILEPATH/livebirths_by_sex.csv")
setnames(birth_forecast, "population", "value")
#create the totals for the combined sex group for the birth forecasts:
aggregate_pop_sex<- function(dt_233){
  dt_233[, pop_both_sex:= sum(value), by = c("location_id", "year_id", "scenario", "age_group_id")]
  
  agg_df_sex<- unique(dt_233[,-c("sex_id", "value")])
  agg_df_sex$sex_id<- 3
  setnames(agg_df_sex, "pop_both_sex", "value")
  
  agg_df<- rbind(dt_233[,-"pop_both_sex"], agg_df_sex)
  
  
  
  return(agg_df)
}
birth_forecast_agg_sex<- aggregate_pop_sex(birth_forecast)

#Open age group 34 from MR-BRT modeling of population:
pop_34_forecast<- readRDS("FILEPATH/pop34_locs20240130.RDS")
pop_34_forecast[, scenario:= 0]

#Open age group 388 from MR-BRT modeling of population:
pop_388_forecast<- readRDS(paste0("FILEPATH/pop388_locs20240130.RDS"))
pop_388_forecast[, scenario:= 0]

#Aggregate populations to the desired the location hierarchy:
aggregate_pop_loc<- function(dt_204){
  
  loc.ref <- get_location_metadata(location_set_id = 92, release_id = 9)
  #all countries
  agg_df <- merge(loc.ref[level==3, c("location_id","location_name","region_id", "super_region_id")], dt_204, by = c("location_id") )
  
  agg_df[, pop_r:= sum(value), by = c("region_id", "year_id", "sex_id", "age_group_id")]
  agg_df[, pop_sr:= sum(value), by = c("super_region_id", "year_id", "sex_id", "age_group_id")]
  agg_df[, pop_g:= sum(value), by = c("year_id", "sex_id", "age_group_id")]
  
  #create region specic data frame
  agg_df_r<- unique(agg_df[,-c("location_id","location_name", "super_region_id","value", "pop_sr", "pop_g")])
  setnames(agg_df_r, c("region_id", "pop_r"), c("location_id", "value"))
  
  
  #create super region specic data frame
  agg_df_sr<- unique(agg_df[,-c("location_id","location_name", "region_id","value", "pop_r", "pop_g")])
  setnames(agg_df_sr, c("super_region_id", "pop_sr"), c("location_id", "value"))
  
  
  #create global data frame
  agg_df_g<- unique(agg_df[,-c("location_id","location_name", "region_id", "super_region_id","value", "pop_sr", "pop_r")])
  agg_df_g$location_id<- 1  
  setnames(agg_df_g, "pop_g", "value")
  
  
  return(rbind(agg_df_g, agg_df_sr, agg_df_r, dt_204))
}
pop_34_forecast_agg<- aggregate_pop_loc(pop_34_forecast[scenario==0 & year_id>= 2021])
pop_388_forecast_agg<- aggregate_pop_loc(pop_388_forecast[scenario==0 & year_id>= 2021])

#intercept shift the population forecasts for age groups 34 and 388:
pop_388_forecast_shifted<- pop_intercept_shift(pop_388_forecast_agg)
pop_34_forecast_shifted<- pop_intercept_shift(pop_34_forecast_agg)

#rbind the population forecast with the sex-combined birth forecasts and shifted+aggregated population forecasts for age groups 34 and 388:
pop_forecast_shifted<- rbind(pop_forecast, birth_forecast_agg_sex[,-"age_group_name"], pop_388_forecast_shifted, pop_34_forecast_shifted)
#Sum the component populations to create the new age group 390 (under 6 months):
pop_forecast_shifted_390<- pop_forecast_shifted[year_id>2021 & age_group_id%in%c(388, 2, 3)]
pop_forecast_shifted_390[, pop_tot:= sum(value), by= c("location_id", "sex_id", "year_id", "scenario")]
pop_forecast_shifted_390<- pop_forecast_shifted_390[age_group_id==388, -"value"]
pop_forecast_shifted_390[, age_group_id:= 390]
setnames(pop_forecast_shifted_390, "pop_tot", "value")
#merge the new age group 390 into the population forecast:
pop_forecast_shifted<- rbind(pop_forecast_shifted, pop_forecast_shifted_390)

#Save the shifted population forecast that includes new age group ids:
saveRDS(pop_forecast_shifted, paste0("FILEPATH",gsub("-", "", Sys.Date()),"_pop_reshifted_34_390.rds"))