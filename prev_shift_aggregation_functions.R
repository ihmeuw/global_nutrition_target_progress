################################################################################
## DESCRIPTION: Shift the forecasted prevalence posterior draws from 2021 onwards by minimizing the logit-space differences between the year 2021 forecast draws and the 2021 GBD estimate.
#Then aggregate the shifted forecast draws to the desired location, sex, age groups, and then calculate the mean prevalence and uncertainty interval for every location-year.
## INPUTS: 
#a) Separate data frames containing forecasted country-age-sex-year prevalence draws from 2021 to 2050 for each of the GNT indicators. 
#b) GBD 2021 country-age-sex specific mean population estimates for year 2021 in age group ids 164, 1, 2,3,4, 5, 388, 8:14, 24, 34
#c) Data frame containing country-age-sex-year specific mean population forecasts from 2021 to 2050 for all age groups relevant to GNT indicators 
## OUTPUTS: 
#a) Separate data frames for each indicator containing forecasted mean prevalence and uncertainty intervals from 2021 to 2050 that have been intercept shifted 
# to the final year of GBD 2021 and aggregated to relevant sex and age groups for the GNTs.

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
#create age group 390
pop_390_2021<- pop_2021[age_group_id%in%c(388, 2, 3)]
pop_390_2021[, pop_tot:= sum(value), by= c("location_id", "sex_id", "year_id")]
pop_390_2021<- pop_390_2021[age_group_id==388, -"value"]
pop_390_2021[, age_group_id:= 390]
setnames(pop_390_2021, "pop_tot", "value")
pop_2021<- rbind(pop_2021, pop_390_2021)

#Open the new shifted population forecasts that includes new age group ids:
pop_forecast_shifted<- readRDS("FILEPATH/20240130_pop_reshifted_34_390.rds")

#Function to format the age-aggregated forecast draws, and calculate mean estimates, lower, and upper. By location, year, sex, and scenario.
prepare_forecast_format <- function(df, standard.locs){
  
  loc.meta <- standard.locs[, c("location_id", "location_name", "location_type")]
  
  df <- merge(df, loc.meta, by = "location_id")
  
  df[sex_id == 1, sex := "Male"]
  df[sex_id == 2, sex := "Female"]
  df[sex_id == 3, sex := "Both"]
  df[, lys := paste0(location_id, "_", year_id, "_", sex_id, "_", scenario)]
  df[, mean_forecast:= mean(pred_prev_cascade), by = lys]
  df[, lower_forecast:= quantile(pred_prev_cascade, probs = .025), by = lys]
  df[, upper_forecast:= quantile(pred_prev_cascade, probs = .975), by = lys]
  df <- df[draw == 0]
  df$draw <- NULL
  df$pred_prev <- NULL
  df$pred_prev_cascade<- NULL
  
  return(df)
}
#Function to format the age-aggregated logit-space forecast draws, and calculate mean estimates, lower, and upper. By location, year, sex, and scenario.
prepare_forecast_format_logit <- function(df, standard.locs){
  
  loc.meta <- standard.locs[, c("location_id", "location_name", "location_type")]
  
  df <- merge(df, loc.meta, by = "location_id")
  
  df[sex_id == 1, sex := "Male"]
  df[sex_id == 2, sex := "Female"]
  df[sex_id == 3, sex := "Both"]
  df[, lys := paste0(location_id, "_", year_id, "_", sex_id, "_", scenario)]
  df[, mean_logit_forecast:= mean(logit(pred_prev_cascade)), by = lys]
  df[, lower_logit_forecast:= quantile(logit(pred_prev_cascade), probs = .025), by = lys]
  df[, upper_logit_forecast:= quantile(logit(pred_prev_cascade), probs = .975), by = lys]
  df <- df[draw == 0]
  df$draw <- NULL
  df$pred_prev <- NULL
  df$pred_prev_cascade<- NULL
  
  return(df)
}
#Function to format the forecasted draws, and calculate mean estimates, lower, and upper. By location, year, age, sex, and scenario.
prepare_df_format <- function(df, standard.locs){
  df[age_group_id == 1, age_group_name := "Under 5"]
  df[age_group_id == 2, age_group_name := "Early Neonatal"]
  df[age_group_id == 3, age_group_name := "Late Neonatal"]
  df[age_group_id == 4, age_group_name := "Post Neonatal"]
  df[age_group_id == 5, age_group_name := "1 to 4"]
  df[age_group_id == 6, age_group_name := "5 to 9"]
  df[age_group_id == 7, age_group_name := "10 to 14"]
  df[age_group_id == 8, age_group_name := "15 to 19"]
  df[age_group_id == 9, age_group_name := "20 to 24"]
  df[age_group_id == 10, age_group_name := "25 to 29"]
  df[age_group_id == 11, age_group_name := "30 to 34"]
  df[age_group_id == 12, age_group_name := "35 to 39"]
  df[age_group_id == 13, age_group_name := "40 to 44"]
  df[age_group_id == 14, age_group_name := "45 to 49"]
  df[age_group_id == 15, age_group_name := "50 to 54"]
  df[age_group_id == 16, age_group_name := "55 to 59"]
  df[age_group_id == 17, age_group_name := "60 to 64"]
  df[age_group_id == 18, age_group_name := "65 to 69"]
  df[age_group_id == 19, age_group_name := "70 to 74"]
  df[age_group_id == 20, age_group_name := "75 to 79"]
  df[age_group_id == 30, age_group_name := "80 to 84"]
  df[age_group_id == 31, age_group_name := "85 to 89"]
  df[age_group_id == 32, age_group_name := "90 to 94"]
  df[age_group_id == 235, age_group_name := "95 plus"]
  df[age_group_id == 257, age_group_name := "20-79 years"]
  df[age_group_id == 388, age_group_name := "1-5 Months"]
  df[age_group_id == 390, age_group_name := "<6 months"]
  df[age_group_id == 389, age_group_name := "6-11 Months"]
  df[age_group_id == 238, age_group_name := "12-23 Months"]
  df[age_group_id == 34, age_group_name := "2-4 Years"]
  df[age_group_id == 24, age_group_name := "15-49 years"]
  df[age_group_id == 164, age_group_name := "Birth"]
  
  df[rei_id== 241, indicator:="stunting"]
  df[rei_id== 240, indicator:="wasting"]
  df[rei_id== 371, indicator:="overweight"]
  df[rei_id== 335, indicator:="LBW"]
  df[rei_id== 192, indicator:="anemia"]
  df[rei_id== 136, indicator:="EBF"]
  #df[covariate_id== 2334, indicator:= "Alcohol (liters)"]
  
  loc.meta <- standard.locs[, c("location_id", "location_name", "location_type")]
  
  df <- merge(df, loc.meta, by = "location_id")
  
  #df[measure_id == 5, measure_name := "Prevalence"]
  #df[metric_id == 3, metric_name := "Rate"]
  df[sex_id == 1, sex := "Male"]
  df[sex_id == 2, sex := "Female"]
  df[sex_id == 3, sex := "Both"]
  #df[, lyas := paste0(location_id, "_", year_id, "_", age_group_id, "_", sex_id)]
  df[, lyass := paste0(location_id, "_", year_id, "_", age_group_id, "_", sex_id, "_", scenario)]
  
  df[, mean := mean(val), by = lyass]
  df[, lower:= quantile(val, probs = .025), by = lyass]
  df[, upper:= quantile(val, probs = .975), by = lyass]
  df <- df[draw == 0]
  df$val <- NULL
  
  return(df)
  
}


intercept_shift<- function(rei, agid, sid, scenario){
forecast<- readRDS(paste0("FILEPATH/", rei,"_", agid, "_", sid, "_", scenario,".RDS"))
#down-sample to the first 500 draws, if applicable:
forecast<- forecast[draw < 500, ]
gbd_estimate_2021 <- readRDS(paste0("FILEPATH/", rei, "_y2021.RDS"))
#determine the difference between the mean forecast for Y2021 and the mean prevalence estimate from GBD 2021
forecast_2021_estimate<- prepare_forecast_format(forecast[year_id==2021,], loc.met)#forecast[year_id==2021,]#prepare_forecast_format(forecast[year_id==2021,], loc.met)
temp<- merge(forecast_2021_estimate, gbd_estimate_2021[age_group_id==agid & sex_id== sid, c("location_id","val")], by=c("location_id"))
temp_draw<- merge(forecast[year_id==2021 & location_id%in%locs,], gbd_estimate_2021[age_group_id==agid & sex_id== sid, c("location_id","val")], by=c("location_id"))

#Calculate intercept shift for each forecast draw, in logit-space to ensure that shifted natural-space prevalence value does not go below zero:
error_function<- function(k){
  y_shift= inv.logit(logit(y) + k)
  return((mean(y_shift)- y_bar_new)^2)
}


for(loc in loc.met[level==3,]$location_id){
y_bar_new<- unique(temp_draw[location_id==loc,]$val)
y<- temp_draw[location_id==loc,]$pred_prev_cascade

fit.k <- optim(1, error_function, method = "L-BFGS-B", lower = -10, upper= 10)
#this is the intercept shift to use for all draws of this location:
if(loc==33){
  temp_draw$shift_logit<- .99999 
  temp_draw[!(location_id%in%loc.met[level==3,]$location_id)]$shift_logit<- NA

}
temp_draw[location_id==loc, shift_logit:= unlist(fit.k$par)]
}

temp[, shift_logit:= logit(val) - logit(mean_forecast)] 

setnames(temp_draw, "val", "est_2021")

#create a new data frame with shift values of the logit-space difference of the mean 2021 estimate and the 2021 draw value 
forecast_shifted<- merge(forecast, temp_draw[draw==1,c("location_id", "shift_logit", "est_2021")], by= "location_id", all.x = TRUE)

#Add the shift in logit space to the logit-space forecast draws, then take the inverse logit to get the shifted natural-space forecast draws
forecast_shifted[, val:= inv.logit(logit(pred_prev_cascade) + shift_logit)] 

if(rei== 136){agid<- 390}
forecast_shifted$age_group_id<- agid

#compare the shifted y2021 mean forecast to the GBD2021 forecast for that year
test<- data.table(prepare_df_format(forecast_shifted, loc.met))
test[, est_diff:= mean - est_2021]
summary(test[year_id==2021]$est_diff)

#merge in 2021+ population estimates
forecast_shifted<- merge(forecast_shifted, pop_forecast_shifted[year_id>=2021 & scenario==0 & age_group_id==agid & sex_id== sid, c("value", "year_id", "location_id")], by= c("location_id", "year_id"), all.x = TRUE)
forecast_shifted[, population:= value]
 return(forecast_shifted[,-c("pred_prev", "pred_prev_cascade", "value", "shift_logit", "est_2021", "age_group_name", "indicator")])#"shift_difference")])

}

#intercept-shift the prevalence forecast draws to final GBD estimates at the finest location, year, age, sex (and scenario) level
lbw<- rbind(intercept_shift(rei = 335, agid = 164, sid = 1, scenario = 0), 
            intercept_shift(rei = 335, agid = 164, sid = 2, scenario = 0))

ebf<- rbind(intercept_shift(rei = 136, agid = 3, sid = 1, scenario = 0), 
            intercept_shift(rei = 136, agid = 3, sid = 2, scenario = 0))

waste<- rbind(intercept_shift(rei = 240, agid = 2, sid = 1, scenario = 0), 
              intercept_shift(rei = 240, agid = 2, sid = 2, scenario = 0), 
              intercept_shift(rei = 240, agid = 3, sid = 1, scenario = 0), 
              intercept_shift(rei = 240, agid = 3, sid = 2, scenario = 0), 
              intercept_shift(rei = 240, agid = 4, sid = 1, scenario = 0),
              intercept_shift(rei = 240, agid = 4, sid = 2, scenario = 0), 
              intercept_shift(rei = 240, agid = 5, sid = 1, scenario = 0), 
              intercept_shift(rei = 240, agid = 5, sid = 2, scenario = 0))

stunt<- rbind(intercept_shift(rei = 241, agid = 2, sid = 1, scenario = 0), 
              intercept_shift(rei = 241, agid = 2, sid = 2, scenario = 0), 
              intercept_shift(rei = 241, agid = 3, sid = 1, scenario = 0), 
              intercept_shift(rei = 241, agid = 3, sid = 2, scenario = 0), 
              intercept_shift(rei = 241, agid = 4, sid = 1, scenario = 0), 
              intercept_shift(rei = 241, agid = 4, sid = 2, scenario = 0), 
              intercept_shift(rei = 241, agid = 5, sid = 1, scenario = 0), 
              intercept_shift(rei = 241, agid = 5, sid = 2, scenario = 0))

overweight<- rbind(intercept_shift(rei = 371, agid = 34, sid = 1, scenario = 0), 
                   intercept_shift(rei = 371, agid = 34, sid = 2, scenario = 0))

anemia<-  rbind(intercept_shift(rei = 95, agid = 8, sid = 2, scenario = 0), 
                intercept_shift(rei = 95, agid = 9, sid = 2, scenario = 0), 
                intercept_shift(rei = 95, agid = 10, sid = 2, scenario = 0), 
                intercept_shift(rei = 95, agid = 11, sid = 2, scenario = 0),
                intercept_shift(rei = 95, agid = 12, sid = 2, scenario = 0), 
                intercept_shift(rei = 95, agid = 13, sid = 2, scenario = 0), 
                intercept_shift(rei = 95, agid = 14, sid = 2, scenario = 0))

#aggregate the prevalence forecast draws up the GBD location hierarchy
aggregate_draws_loc<- function(shifted_forecast_df){

loc.ref <- get_location_metadata(location_set_id = 92, release_id = 9)
#all countries
agg_df <- merge(loc.ref[level==3, c("location_id","location_name","region_id", "super_region_id")], shifted_forecast_df[year_id>2020,], by = c("location_id") )

agg_df[, wt_r_val:= weighted.mean(val, population), by = c("region_id", "year_id", "sex_id", "age_group_id", "scenario", "draw")]
agg_df[, wt_sr_val:= weighted.mean(val, population), by = c("super_region_id", "year_id", "sex_id", "age_group_id", "scenario", "draw")]
agg_df[, wt_g_val:= weighted.mean(val, population), by = c("year_id", "sex_id", "age_group_id", "scenario", "draw")]

#create region specic data frame
agg_df_r<- unique(agg_df[,-c("location_id","location_name", "super_region_id","population", "val", "wt_sr_val", "wt_g_val")])
setnames(agg_df_r, "region_id", "location_id")
agg_df_r<- merge(agg_df_r, pop_2021[,c("age_group_id", "sex_id", "population", "year_id", "location_id")], by= c("location_id", "year_id", "age_group_id", "sex_id"), all.x = TRUE)
agg_df_r<- merge(agg_df_r, pop_forecast_shifted[year_id>2021 & scenario==0, c("location_id", "year_id", "age_group_id", "sex_id", "value")], by= c("location_id", "year_id", "age_group_id", "sex_id"), all.x = TRUE)
agg_df_r[year_id>2021, population:= value]
agg_df_r<- agg_df_r[,-c("value")]
setnames(agg_df_r, "wt_r_val", "val")

#create super region specic data frame
agg_df_sr<- unique(agg_df[,-c("location_id","location_name", "region_id","population", "val", "wt_r_val", "wt_g_val")])
setnames(agg_df_sr, "super_region_id", "location_id")
agg_df_sr<- merge(agg_df_sr, pop_2021[,c("age_group_id", "sex_id", "population", "year_id", "location_id")], by= c("location_id", "year_id", "age_group_id", "sex_id"), all.x = TRUE)
agg_df_sr<- merge(agg_df_sr, pop_forecast_shifted[year_id>2021 & scenario==0, c("location_id", "year_id", "age_group_id", "sex_id", "value")], by= c("location_id", "year_id", "age_group_id", "sex_id"), all.x = TRUE)
agg_df_sr[year_id>2021, population:= value]
agg_df_sr<- agg_df_sr[,-c("value")]
setnames(agg_df_sr, "wt_sr_val", "val")

#create global data frame
agg_df_g<- unique(agg_df[,-c("location_id","location_name", "region_id", "super_region_id", "population", "val", "wt_r_val", "wt_sr_val")])
agg_df_g[, location_id:= 1]
agg_df_g<- merge(agg_df_g, pop_2021[,c("age_group_id", "sex_id", "population", "year_id", "location_id")], by= c("location_id", "year_id", "age_group_id", "sex_id"), all.x = TRUE)
agg_df_g<- merge(agg_df_g, pop_forecast_shifted[year_id>2021 & scenario==0, c("location_id", "year_id", "age_group_id", "sex_id", "value")], by= c("location_id", "year_id", "age_group_id", "sex_id"), all.x = TRUE)
agg_df_g[year_id>2021, population:= value]
agg_df_g<- agg_df_g[,-c("value")]
setnames(agg_df_g, "wt_g_val", "val")

return(rbind(agg_df_g, agg_df_sr, agg_df_r, shifted_forecast_df))
}

lbw_loc_agg<- aggregate_draws_loc(lbw)
ebf_loc_agg<- aggregate_draws_loc(ebf)
waste_loc_agg<- aggregate_draws_loc(waste)
stunt_loc_agg<- aggregate_draws_loc(stunt)
overweight_loc_agg<- aggregate_draws_loc(overweight)
anemia_loc_agg<- aggregate_draws_loc(anemia)

#aggregate the prevalence forecast draws to combine both-sexes
aggregate_draws_sex<- function(shifted_forecast_233_df){

  shifted_forecast_233_df[, wt_sex_val:= weighted.mean(val, population), by = c("location_id", "year_id", "scenario", "age_group_id", "draw")]
  
  agg_df_sex<- unique(shifted_forecast_233_df[,-c("sex_id", "val", "population")])
  agg_df_sex$sex_id<- 3
  agg_df_sex<- merge(agg_df_sex, pop_2021[,c("age_group_id", "sex_id", "population", "year_id", "location_id")], by= c("location_id", "year_id", "age_group_id", "sex_id"), all.x = TRUE)
  agg_df_sex<- merge(agg_df_sex, pop_forecast_shifted[year_id>2021 & scenario==0, c("location_id", "year_id", "age_group_id", "sex_id", "value")], by= c("location_id", "year_id", "age_group_id", "sex_id"), all.x = TRUE)
  agg_df_sex[year_id>2021, population:= value]
  agg_df_sex<- agg_df_sex[,-c("value")]
  setnames(agg_df_sex, "wt_sex_val", "val")
  
  agg_df<- rbind(shifted_forecast_233_df[,-"wt_sex_val"], agg_df_sex)
  
  return(agg_df)
  
}

lbw_sex_agg<- aggregate_draws_sex(lbw_loc_agg)
ebf_sex_agg<- aggregate_draws_sex(ebf_loc_agg)
waste_sex_agg<- aggregate_draws_sex(waste_loc_agg)
stunt_sex_agg<- aggregate_draws_sex(stunt_loc_agg)
overweight_sex_agg<- aggregate_draws_sex(overweight_loc_agg)

#aggregate the prevalence forecast draws to desired age groups
aggregate_draws_age<- function(sex_agg_forecast_df, agid_out){
  sex_agg_forecast_df[, wt_age_val:= weighted.mean(val, population), by = c("location_id", "year_id", "scenario", "sex_id", "draw")]
  
  agg_df_age<- unique(sex_agg_forecast_df[,-c("age_group_id", "val", "population")])
  agg_df_age$age_group_id<- agid_out
  agg_df_age<- merge(agg_df_age, pop_2021[,c("age_group_id", "sex_id", "population", "year_id", "location_id")], by= c("location_id", "year_id", "age_group_id", "sex_id"), all.x = TRUE)
  agg_df_age<- merge(agg_df_age, pop_forecast_shifted[year_id>2021 & scenario==0, c("location_id", "year_id", "age_group_id", "sex_id", "value")], by= c("location_id", "year_id", "age_group_id", "sex_id"), all.x = TRUE)
  agg_df_age[year_id>2021, population:= value]
  agg_df_age<- agg_df_age[,-c("value")]
  setnames(agg_df_age, "wt_age_val", "val")
  
  agg_df<- rbind(sex_agg_forecast_df[, -"wt_age_val"], agg_df_age)
  
  return(agg_df)
  
}

lbw_age_agg<- aggregate_draws_age(lbw_sex_agg, 164)
waste_age_agg<- aggregate_draws_age(waste_sex_agg, 1)
stunt_age_agg<- aggregate_draws_age(stunt_sex_agg, 1)
anemia_age_agg<- aggregate_draws_age(anemia_loc_agg, 24) #anemia did not need sex aggregation because females-only

#Calculate the mean prevalence and uncertainty interval for each of the forecasted indicators
lbw_prev_forecast<- prepare_df_format(lbw_sex_agg, loc.met)
ebf_prev_forecast<- prepare_df_format(ebf_sex_agg, loc.met) #EBF did not need age-group aggregation
waste_prev_forecast<- prepare_df_format(waste_age_agg, loc.met)
stunt_prev_forecast<- prepare_df_format(stunt_age_agg, loc.met)
overweight_prev_forecast<- prepare_df_format(overweight_sex_agg, loc.met) #overweight did not need age-group aggregation, because already the desired age group
anemia_prev_forecast<- prepare_df_format(anemia_age_agg, loc.met)

#Save the forecasted prevalence data frames
saveRDS(lbw_prev_forecast,paste0("FILEPATH/lbw_forecast.RDS"))
saveRDS(ebf_prev_forecast,paste0("FILEPATH/ebf_forecast.RDS"))
saveRDS(stunt_prev_forecast,paste0("FILEPATH/stunt_forecast.RDS"))
saveRDS(waste_prev_forecast,paste0("FILEPATH/waste_forecast.RDS"))
saveRDS(overweight_prev_forecast,paste0("FILEPATH/overweight_forecast.RDS"))
saveRDS(anemia_prev_forecast,paste0("FILEPATH/anemia_forecast.RDS"))

