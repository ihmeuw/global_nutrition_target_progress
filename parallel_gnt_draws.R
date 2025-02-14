################################################################################
## DESCRIPTION: Parallel script for array job to run MR-BRT and cascading splines models to predict prevalence of Global Nutrition Target (GNT) indicators from Summary Exposure Values and socio-demographic index
## INPUTS: 
#a) Forecasted country-age-sex-year Summary Exposure Value draws for each GNT indicator or related impairment (i.e., iron deficiency, in the case of anaemia) 
#b) Forecasts of country-level socio-demographic index
## OUTPUTS: 
#a) Predicted country-age-sex-year prevalence draws from 2021 to 2050 

#Before running this script, follow instructions here: https://github.com/ihmeuw-msca/mrtoolR
##############################################################################################################################################
##############################################################################################################################################
########################################################### SET UP ###########################################################################
##############################################################################################################################################
##############################################################################################################################################
library(reticulate)
reticulate::use_python("FILEPATH")
library(mrbrt001, lib.loc = "FILEPATH")
mr <- import("mrtool")
library(dplyr)
library(data.table)
library(boot)

invisible(sapply(list.files("/share/cc_resources/libraries/current/r/", full.names = T), source))


#pull information from param_map
if(interactive()){
  task_id<- 
  param_map_filepath<- "FILEPATH/param_map.csv"
}else{
  args <- commandArgs(trailingOnly = TRUE) # you're telling this script that there is important information at the end of the qsub you launched
  param_map_filepath <- args[1] # you attached the param map filepath to the end of the qsub
  task_id <- ifelse(Sys.getenv('SLURM_ARRAY_TASK_ID') != '', as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID')), NA)# this script knows which row of the param map it corresponds to
}
param_map <- fread(param_map_filepath)

print(task_id)
print(param_map[task_id])

r = param_map[task_id, rei_id] 
agid = param_map[task_id, age_group_id] 
sid = param_map[task_id, sex_id] 
scenid = param_map[task_id, scenario] 

output_dir <- file.path("FILEPATH")
output_dir2<-"FILEPATH"
if (r==335 | r== 136 | r== 241 | r== 240 | r== 95){
   sev_df <- readRDS(file.path("FILEPATH/SEV_prev2021.RDS"))#retrospective SEV and prevalence estimates
   
  if (r== 335) {
    ind<-"lbw"
    mod_forecast <- readRDS("FILEPATH/LBW_sev_forecast.RDS")
    standard <- "mod_full"
    cascade <- "mod_cascade_full"
    #change the age_group_id to match LBW's 164
    sev_df[rei_id==335, age_group_id:=164]
  }
  if (r== 136) {
    #agid <- 3
    ind<- "ebf"
    mod_forecast <- readRDS("FILEPATH/EBF_sev_forecast.RDS")
    mod_forecast <- mod_forecast[,-"location_name"]
    standard <- "mod_full"
    cascade <- "mod_cascade_full"
  }
  if (r== 241) {
    ind<- "stunt"
    mod_forecast <- readRDS("FILEPATH/stunt_sev_forecast.RDS")
    mod_forecast <- mod_forecast[,-"location_name"]
    standard <- "mod_full"
    cascade <- "mod_cascade_full"
  }
  if (r== 240) {
    ind<- "waste"
    mod_forecast <- readRDS("FILEPATH/waste_sev_forecast.RDS")
    mod_forecast <- mod_forecast[,-"location_name"]
    standard <- "mod_full"
    cascade <- "mod_cascade_full"
  }
  if (r== 95) {
    ind<- "anem"
    sev_df <- readRDS(file.path("FILEPATH/anem_prev2021.RDS"))
    mod_forecast <- readRDS(paste0(output_dir,"FILEPATH/anem_sev_forecast.RDS"))
    standard <- "mod_full"
    cascade <- "mod_cascade_full"
  }

  fit_global_full <- function(model_data, r, agid, sid) {
    #subset to only contain both sexes and split into training and testing dataframes
    train_df<- model_data[rei_id==r & sex_id== sid & age_group_id== agid,]

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
                       spline_knots = array(seq(0, 1, by = 0.25)), # put a spline every 0.25
                       spline_degree = 3L,
                       spline_knots_type = 'frequency',
                       spline_r_linear = TRUE,
                       spline_l_linear = FALSE
        )
      ))
    stand_mod$fit_model(inner_print_level = 5L, inner_max_iter = 2000L,  outer_max_iter= 1000L)
    estimated_beta <- stand_mod$summary()[[1]][1, "sdi"]
    stand_mod$cov_models[[which(stand_mod$cov_names == "sdi")]]$prior_beta_uniform <- matrix(rep(estimated_beta, 2), ncol = 1)
    

    return(stand_mod)
    
  }

    mod_full <- fit_global_full(sev_df,r, agid, sid)
    #pull results from cascading spline from file:
    sex <- ifelse(sid==1, "Male", "Female")
    assign(paste0("mod_cascade_full"), list(working_dir = paste0(output_dir2,  "cascade_prevalence", r, "_", agid, "_", sex)))
}


    if(r!= 371){
      forecast<- as.data.frame(matrix(data=0,nrow = length(unique(mod_forecast$year_id)) * length(unique(mod_forecast$location_id))*500, ncol=6)) #500 draws
    colnames(forecast) <- c("location_id", "year_id","scenario", "draw", "pred_prev", "pred_prev_cascade")    
      
      for(drawid in c(0:499)) { 
      
        forecast_df <- mod_forecast[rei_id== r & sex_id== sid & age_group_id== agid & scenario == scenid & draw== drawid,]

        #prediction data- forecast
        dat_pred_forecast <- MRData()
        dat_pred_forecast$load_df(
          data = forecast_df,
          col_covs=list("sdi", "sev"), col_study_id = "location_id")
        
        #Forecasting the future
        pred_logit <- get(standard)$predict(data = dat_pred_forecast, predict_for_study= TRUE, sort_by_data_id= TRUE)
        preds_loc_forecast <- forecast_df
        preds_loc_forecast$pred_logit <- pred_logit
        preds_loc_forecast[,pred_prev:= inv.logit(pred_logit)]
        
        #forecasting the future
        preds_loc_forecast <- as.data.table(predict_spline_cascade(fit = get(cascade), newdata = preds_loc_forecast))
        preds_loc_forecast[, pred_prev_cascade := inv.logit(pred)]
        
        #add scenario
        preds_loc_forecast[,scenario:=scenid]
        #add draw
        preds_loc_forecast[,draw:=drawid]
        
        #add the relevant row numbers
        forecast[((length(unique(mod_forecast$year_id)) * length(unique(mod_forecast$location_id))*drawid)+1):((drawid+1) * length(unique(mod_forecast$year_id)) * length(unique(mod_forecast$location_id))), ]<- preds_loc_forecast[,c("location_id", "year_id","scenario", "draw", "pred_prev", "pred_prev_cascade")]
      }
    }

    
    #Childhood overweight     
    if(r== 371){
      
        overweight_df <- readRDS(file.path("FILEPATH/overweight_prev2021.RDS"))
        overweight_df<- overweight_df[year_id<2021,]
        mod_forecast <- readRDS("FILEPATH/overweight_sev_forecast.RDS")
        standard <- "mod_full"
        cascade <- "mod_cascade_full"

        agid <- 34
        
        fit_global_full <- function(model_data, r, agid, sid) {
          #subset to only contain both sexes and split into training and testing dataframes
          train_df<- model_data[rei_id==r & sex_id== sid & age_group_id== agid,]
          
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
                             spline_knots = array(seq(0, 1, by = 0.25)), # put a spline every 0.25
                             spline_degree = 3L,
                             spline_knots_type = 'frequency',
                             spline_r_linear = TRUE,
                             spline_l_linear = FALSE#,
              )
            ))
          
          stand_mod$fit_model(inner_print_level = 5L, inner_max_iter = 1000L,  outer_max_iter= 500L)

          
          return(stand_mod)
          
        }
        
        mod_full <- fit_global_full(overweight_df,r, agid, sid)
        
        #pre-run the overweight cascade prevalence models
        sex <- ifelse(sid==1, "Male", "Female")
        assign(paste0("mod_cascade_full"), list(working_dir = paste0(output_dir2,  "cascade_prevalence", r, "_", agid, "_", sex)))
        
  forecast<- as.data.frame(matrix(data=0,nrow = length(unique(mod_forecast$year_id)) * length(unique(mod_forecast$location_id))*500, ncol=6))#nrow = 6324*1000, ncol=6))
  colnames(forecast) <- c("location_id", "year_id","scenario", "draw", "pred_prev", "pred_prev_cascade")
      
    for(drawid in c(0:499)) {
 
      forecast_df <- mod_forecast[rei_id== r & sex_id== sid & age_group_id== agid & scenario == scenid & draw== drawid,]
      
      #prediction data- forecast
      dat_pred_forecast <- MRData()
      dat_pred_forecast$load_df(
        data = forecast_df,
        col_covs=list("sdi", "sev_kids"), col_study_id = "location_id")
      
      #Forecasting the future
      pred_logit <- get(standard)$predict(data = dat_pred_forecast, predict_for_study= TRUE, sort_by_data_id= TRUE)
      preds_loc_forecast <- forecast_df
      preds_loc_forecast$pred_logit <- pred_logit
      preds_loc_forecast[,pred_prev:= inv.logit(pred_logit)]
      
      #forecasting the future
      preds_loc_forecast <- as.data.table(predict_spline_cascade(fit = get(cascade), newdata = preds_loc_forecast))
      preds_loc_forecast[, pred_prev_cascade := inv.logit(pred)]
      
      #add scenario
      preds_loc_forecast[,scenario:=scenid]
      #add draw
      preds_loc_forecast[,draw:=drawid]

      forecast[((length(unique(mod_forecast$year_id)) * length(unique(mod_forecast$location_id))*drawid)+1):((drawid+1) * length(unique(mod_forecast$year_id)) * length(unique(mod_forecast$location_id))), ]<- preds_loc_forecast[,c("location_id", "year_id","scenario", "draw", "pred_prev", "pred_prev_cascade")]
      
    }
    }
    forecast<-as.data.table(forecast)
    forecast[, sex_id:=sid]
    forecast[,rei_id:= r]
    forecast[,location_id:= as.numeric(location_id)]

saveRDS(forecast, paste0("FILEPATH/", r, "_", agid, "_", sid, "_", scenid, ".RDS"))