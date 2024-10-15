################################################################################
## DESCRIPTION: Use meta-regression Bayesian prior tool (MR-BRT) and cascading splines to model population of children 1-5 months of age (age_group_id 388)
## INPUTS: 
  #a) Separate population estimates of children in Late Neonatal (age_group_id 3), post-neonatal period (age_group_id 4), and country-level socio-demographic index (SDI) 
  #b) Separate population forecasts of children in Late Neonatal (age_group_id 3), post-neonatal period (age_group_id 4), and country-level SDI
## OUTPUTS: 
  #a) Predicted population of age 1-5 months in training sample (1990-2014), predicted predicted population in testing sample (2015-2021), root-mean squared error (RMSE) ##
  #b) Predicted population forecasts of age 1-5 months from 2021 to 2050 ##
## AUTHOR: 
## DATE:
################################################################################

#predicting population of children 1-5 months of age
library(tidyverse)
library(haven)
library(parallel)
library(tidyselect)
invisible(sapply(list.files("/share/cc_resources/libraries/current/r/", full.names = T), source))

invisible(sapply(list.files("FILEPATH", full.names = T), source))

path <- paste0("FILEPATH")

#MR-BRT and cascading splines functions require docker. Instructions are available from: https://github.com/ihmeuw-msca/mrtoolr/tree/main
library(reticulate); use_condaenv("mrtool-0.0.1"); library(mrtoolr)
mr <- import("mrtoolr")

# Directories -------------------------------------------------------------
output_dir <- file.path("FILEPATH")
output_dir1<-"FILEPATH"
package_lib <- sprintf(path)
.libPaths(path)

library(matrixStats)
library(data.table)
library(zoo)
library(plyr)
library(dplyr)
library(scales)
library(ggplot2)

packages <- c("magrittr","ggplot2", "Metrics", "boot", "grid", "gridExtra", "ggpubr", "patchwork")
for(p in packages){
  if(p %in% rownames(installed.packages())==FALSE){
    install.packages(p)
  }
  library(p, character.only = T)
}



loc.met <- get_location_metadata(location_set_id = 92, release_id = 9)
#all countries
loc.met <- loc.met[level==3]
locs <- unique(loc.met$location_id)

locs <- unique(loc.met$location_id)
#Function for formatting data as desired
prepare_df_format <- function(df, standard.locs){
  
  df[age_group_id == 1, age_group_name := "Under 5"]
  df[age_group_id == 2, age_group_name := "Early Neonatal"]
  df[age_group_id == 3, age_group_name := "Late Neonatal"]
  df[age_group_id == 4, age_group_name := "Post Neonatal"]
  
  #for GBD2019 we need to use agegroupid 5 for overweight in chidren
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
  
  df[modelable_entity_id==10556, indicator:="stunting"]
  df[modelable_entity_id==10558, indicator:="wasting"]
  df[modelable_entity_id==20018, indicator:="overweight"]
  df[modelable_entity_id==16282, indicator:="LBW"]
  df[modelable_entity_id==10507, indicator:="anemia"]
  df[modelable_entity_id==20417, indicator:="EBF"]
  df[modelable_entity_id==11231, indicator:="PM2.5"]
  #df[covariate_id== 2334, indicator:= "Alcohol (liters)"]
  
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
id.vars <- c("metric_id", "age_group_id", "location_id", "measure_id", "modelable_entity_id", "sex_id", "year_id", "model_version_id")
hierarchy <- get_location_metadata(92, release_id = 9)

#sdi
sdi <- fread("sdi_retrospective.csv")
setnames(sdi, old=c("mean","lower", "upper"), new=c("sdi","lower_sdi", "upper_sdi"))
#retrospective population estimates 
population <- fread("population_retrospective.csv")
#get the combined population for agegroup 4. This will be used as predictor variable in the cascading spline model
population[age_group_id==3, pop3:= population]
population[, pop3:= na.locf(pop3, na.rm=FALSE), by= c("location_id", "sex_id", "year_id")]
population[age_group_id>3, pop4:= sum(population), by= c("location_id", "sex_id", "year_id")]
population<- population[age_group_id==388,]
population <- merge(population, hierarchy[, .(location_id, super_region_name, region_name)], all.x=T, by = "location_id")
population<- merge(population, sdi, by=c("location_id", "year_id"))

#se based on UI
# population cannot be put into in logit space
population[, se:= (population - lower)/1.96]

saveRDS(population, file.path(output_dir, "ag388_population2021.RDS"))
pop_388_df <- readRDS(file.path(output_dir, "ag388_population2021.RDS"))
sid<- 1
cor(population[sex_id==sid]$pop3, population[sex_id==sid]$population)

cor(population[sex_id==sid]$pop4, population[sex_id==sid]$population)

fit_pop <- function(model_data,sid) {
  #subset to only contain both sexes and split into training and testing dataframes
  train_df<- model_data[year_id>2004 & year_id<2015 & sex_id==sid,]
  test_df <- model_data[year_id>=2015 & sex_id==sid]
  
  #training data
  dat_loc <- MRData()
  dat_loc$load_df(
    data = train_df,
    col_obs = "population", col_obs_se = "se",
    col_covs = list("sdi", "pop3", "pop4"), col_study_id = "location_id"
  )
  #prediction data
  dat_pred <- MRData()
  dat_pred$load_df(
    data = test_df,
    col_covs=list("sdi", "pop3","pop4"), col_study_id = "location_id"
  )
  
  #"standard MRBRT model
  stand_mod <- MRBRT(
    data = dat_loc,
    cov_models = list(
      LinearCovModel("intercept", use_re = TRUE),
      LinearCovModel("sdi", use_re = FALSE),
      LinearCovModel("pop3", use_re = TRUE),
      LinearCovModel("pop4", use_re = FALSE, 
                     use_spline = TRUE,
                     spline_knots = array(seq(0, 1, by = 0.25)), # this tells is to put a spline every 0.25
                     spline_degree = 1L,
                     spline_knots_type = 'frequency',# this specifies to put the knots literally along the domain. Other option is frequency for it to be data density driven
                     spline_r_linear = TRUE,
                     spline_l_linear = FALSE
      )
    ))
  
  stand_mod$fit_model(inner_print_level = 0L, inner_max_iter = 1000L,  outer_max_iter= 500L)
  estimated_beta <- stand_mod$summary()[[1]][1, "sdi"]
  stand_mod$cov_models[[which(stand_mod$cov_names == "sdi")]]$prior_beta_uniform <- matrix(rep(estimated_beta, 2), ncol = 1)
  
  pred_pop <- stand_mod$predict(data = dat_pred, predict_for_study= TRUE, sort_by_data_id= TRUE)
  test_df$pred_pop <- pred_pop
  print(paste0("Population RMSE= ", rmse(test_df$population, test_df$pred_pop)))
  
  return(stand_mod)
  
}

fit_pop_cascade <- function(mod_global, model_data, output_dir, sid) {
  
  train_df<- model_data[year_id>2004 & year_id<2015 & sex_id==sid,]

  model_label_tmp <- paste0("cascade_population", unique(model_data$age_group_id), "_", sid)
  thetas <- c(2,7)
  cascade_fit <- run_spline_cascade(
    stage1_model_object = mod_global, 
    df = train_df, 
    col_obs = "population",
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
  #parLapply(cl, c(1,2), function(sid){
pop_standard<- fit_pop(population, sid)
pop_cascade<- fit_pop_cascade(pop_standard, population, output_dir1, sid)
}






#population dx plots:

custom.col.sr <- c("Central Europe, Eastern Europe, and Central Asia" = "#771155",
                   "High-income" = "#117777",
                   "Latin America and Caribbean" = "#771122",
                   "North Africa and Middle East" = "#777711", 
                   "South Asia" = "#000080",
                   "Southeast Asia, East Asia, and Oceania" = "#114477", 
                   "Sub-Saharan Africa" = "#774411")
custom.col.r <- c("Central Asia" = "#771155",
                  "Central Europe" = "#AA4488",
                  "Eastern Europe" = "#CC99BB",
                  "Australasia" = "#117777",
                  "High-income Asia Pacific" = "#44AAAA",
                  "High-income North America" = "#77CCCC",
                  "Southern Latin America" = "#117744",
                  "Western Europe" = "#44AA77",
                  "Andean Latin America" = "#771122", # red
                  "Caribbean" = "#993344",
                  "Central Latin America" = "#BB5566",
                  "Tropical Latin America" = "#DD7788",
                  "North Africa and Middle East" = "#777711", 
                  "South Asia" = "#000080", # ?
                  "East Asia" = "#114477", 
                  "Oceania" = "#4477AA",
                  "Southeast Asia" = "#77AADD",
                  "Central Sub-Saharan Africa" = "#774411", #orange
                  "Eastern Sub-Saharan Africa" = "#996633",
                  "Southern Sub-Saharan Africa" = "#BB8855",
                  "Western Sub-Saharan Africa" = "#DDAA77")
pdf(file=paste0("FILEPATH", "/mrbrt_plots_pop388_theta2_7_", gsub("-", "", Sys.Date()) ,".pdf"), height = 6.5, width = 10.2)
for(sid in c(1,2)){
  pop_standard<- fit_pop(population, sid)
  pop_cascade<- list(working_dir =  file.path(output_dir1,paste0("cascade_population388", "_", sid)))
  standard <- "pop_standard"
  cascade <- "pop_cascade"
  measure <- "population"
  time_period_in <- "2005 - 2014"
  time_period_out <- "2015 - 2021"
  space <- 0.005
  predictor <- "Age group 4 population"
  exposure <- "population"
  sex <- ifelse(sid==1, "Male", "Female")
  #predictor <- "Iron SEV"
  #exposure <- "sev"
  train_df <- population[year_id>2004 & year_id<2015 & sex_id==sid,]
  test_df <- population[year_id>=2015 & sex_id==sid,]
  #prediction data- in sample
  dat_pred_in <- MRData()
  dat_pred_in$load_df(
    data = train_df,
    col_covs=list("sdi", "pop3","pop4"), col_study_id = "location_id")
  #prediction data- out of sample
  dat_pred_out <- MRData()
  dat_pred_out$load_df(
    data = test_df,
    col_covs=list("sdi", "pop3","pop4"), col_study_id = "location_id")
  #create DF with SEV knots and SDI values at 0.25? intervals
  pop_vec <- seq(min(get(standard)$cov_models[[4]]$spline_knots),max(get(standard)$cov_models[[4]]$spline_knots),by = 5000000)
  num_pop <- length(pop_vec)
  sdi_vec <- array(quantile(population[year_id>2004 & year_id<2015 & sex_id== sid,]$sdi))[2:5]
  knot_df <- data.frame(pop4 = rep(pop_vec, 4), sdi_group = rep(1:4, each = num_pop),
                        sdi = rep(sdi_vec, each = num_pop),
                        location_id = 0, sex_id = sid, pop3 = 3000, se = 1)
  #prediction data- spline plots
  dat_pred_spline <- MRData()
  dat_pred_spline$load_df(
    data = knot_df,
    #col_covs=list("sdi", "sev"), col_study_id = "location_id")
    col_covs=list("sdi", "pop3", "pop4"), col_study_id = "location_id")


pred_spline_pop <- get(standard)$predict(data = dat_pred_spline, predict_for_study= TRUE, sort_by_data_id= TRUE)
preds_spline <- as.data.table(knot_df)
preds_spline$pred_pop <- pred_spline_pop

#plot SEV at certain value of SDI. Y-axis is prevalence (x-value is SEV)
#predict using standard spline model
#In-sample predictive validity
pred_pop <- get(standard)$predict(data = dat_pred_in, predict_for_study= TRUE, sort_by_data_id= TRUE)
preds_loc_in <- train_df
preds_loc_in$pred_pop <- pred_pop

#add residuals
preds_loc_in[, residual:=population- pred_pop ]
preds_loc_in[which.max(preds_loc_in$residual)]$location_name
print(paste0(sex, " population, standard MRBRT (IS) RMSE= ", rmse(preds_loc_in$population, preds_loc_in$pred_pop)))
#Out of sample predictive validity
pred_pop <- get(standard)$predict(data = dat_pred_out, predict_for_study= TRUE, sort_by_data_id= TRUE)
preds_loc_out <- test_df
preds_loc_out$pred_pop <- pred_pop
# df_mod_log <- cbind(get(standard)$data$to_df(), data.frame(w = get(standard)$w_soln))
#add residuals
preds_loc_out[, residual:=population- pred_pop ]
preds_loc_out[which.max(preds_loc_out$residual)]$location_name
print(paste0(sex, " population, standard MRBRT (OOS) RMSE= ", rmse(preds_loc_out$population, preds_loc_out$pred_pop)))
#predict using cascading spline model
#in-sample
preds_loc_in <- as.data.table(predict_spline_cascade(fit = get(cascade), newdata = preds_loc_in))
preds_loc_in[, pred_pop_cascade := pred]
#add residuals
preds_loc_in[, residual_cascade:=population- pred_pop_cascade ]
preds_loc_in[which.max(preds_loc_in$residual_cascade)]$location_name
print(paste0(sex, " population, cascading spline (IS) RMSE= ", rmse(preds_loc_in$population, preds_loc_in$pred_pop_cascade)))
#out of sample
preds_loc_out <- as.data.table(predict_spline_cascade(fit = get(cascade), newdata = preds_loc_out))
preds_loc_out[, pred_pop_cascade := pred]
#add residuals
preds_loc_out[, residual_cascade:=population- pred_pop_cascade ]
preds_loc_out[which.max(preds_loc_out$residual_cascade)]$location_name
print(paste0(sex, " population, cascading spline (OOS) RMSE= ", rmse(preds_loc_out$population, preds_loc_out$pred_pop_cascade)))

# #plot predicted vs. observed population
for(sample in c("IS", "OOS")){
  if(sample=="IS"){
    preds_loc <- preds_loc_in
    time_period <- time_period_in
    space_res <- -1
  }
  if(sample=="OOS"){
    preds_loc <- preds_loc_out
    time_period <- time_period_out
    space_res <- -0.5
  }
 
    title_stand_spline <- paste0(sex, " population: pop in 1-11 months \n vs. predicted  1-5 months ", measure, " : standard MRBRT")
    title_casc_spline <- paste0(sex, " population: pop in 1-11 months \n vs. predicted 1-5 months ", measure, " : cascading spline")
    title_obs_spline <- paste0(sex, " population: pop in 1-11 months \n vs. predicted 1-5 months ", measure, " ")
    title_label_stand <- paste0(sex, " population in 1-5 months predicted vs. observed ", measure, " in ", time_period, ": standard MRBRT")
    title_label_casc <- paste0(sex, " population in 1-5 months predicted vs. observed ", measure, " in ", time_period, ": cascading spline")
    title_label_res_stand<-paste0(sex, " population in 1-5 months residuals over time from ", time_period, ": standard MRBRT")
    title_label_res_casc <- paste0(sex, " population in 1-5 months residuals over time from ", time_period, ": cascading spline")

  gg_stand <- ggplot(preds_loc, aes(y=population, x=pred_pop))+
    geom_line(data = preds_loc,
              aes(color= super_region_name, group = location_id), alpha = 0.3)+
    theme_bw() +
    scale_color_manual(values = custom.col.sr)+
    labs(title = title_label_stand,
         subtitle = paste0("RMSE (", sample,")= ", formatC(rmse(preds_loc$population, preds_loc$pred_pop), digits = 4),
                           ". SDI coefficient = ", formatC(get(standard)$beta_soln[2], digits=3)),
         y = paste("Observed ", measure),
         x = paste("Predicted ", measure),
         color = "Super Region") +
    theme(legend.position="bottom",
          plot.title = element_text(size=8),
          plot.subtitle = element_text(size=7),
          legend.title.align = .5,
          legend.background = element_rect(fill="grey90", size=0.5, linetype="solid", colour ="grey30"),
          legend.key = element_rect(fill = "grey90")) +
    guides(color = guide_legend(override.aes = list(alpha = .5), title.position = "top")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(ylim=c(0, max(preds_loc$population, preds_loc$pred_pop, preds_loc$pred_pop_cascade)), xlim=c(0, max(preds_loc$population, preds_loc$pred_pop, preds_loc$pred_pop_cascade)))+
    #identify the largest 5 residuals from the final year predicted
    geom_text(data=subset(preds_loc, abs(residual)%in%tail(sort( preds_loc[, max(abs(residual)), by = "location_id"]$V1 ),5)) , aes(y=population,x=pred_pop,label=paste0(location_name, " ", formatC(residual, digits=3, format="f") )), nudge_y = space,  size=1.5)
  
  
  
  gg_stand_spline_knots <- ggplot(preds_spline, aes(y=pred_pop, x=get(exposure)))+
    geom_line(data = preds_spline,
              aes(color= sdi, group = sdi_group), alpha = 1, size=1)+
    theme_bw() +
    geom_line(data = preds_loc, aes(y=pred_pop, x=get(exposure), color= sdi, group = location_id), alpha = 0.7)+
    labs(title = title_stand_spline,
         subtitle = paste0("By SDI quartile. RMSE (", sample,")= ", formatC(rmse(preds_loc$population, preds_loc$pred_pop), digits = 4),
                           ". SDI coefficient = ", formatC(get(standard)$beta_soln[2], digits=3)),
         y = paste("Predicted ", measure),
         x = predictor,
         color = "SDI") +
    theme(legend.position="bottom",
          plot.title = element_text(size=8),
          plot.subtitle = element_text(size=7),
          legend.title.align = .5,
          legend.background = element_rect(fill="grey90", size=0.5, linetype="solid", colour ="grey30"),
          legend.key = element_rect(fill = "grey90")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    geom_vline(xintercept = c(get(standard)$cov_models[[3]]$spline_knots) ) +
    coord_cartesian(ylim=c(0, max(preds_loc$population, preds_loc$pred_pop, preds_loc$pred_pop_cascade)))
  
  
  gg_casc_spline_knots <- ggplot(preds_spline, aes(y=pred_pop, x=get(exposure)))+
    geom_line(data = preds_spline,
              aes(color= sdi, group = sdi_group), alpha = 1, size=1)+
    theme_bw() +
    geom_line(data = preds_loc, aes(y=pred_pop_cascade, x=get(exposure), color= sdi, group = location_id), alpha = 0.7)+
    labs(title = title_casc_spline,
         subtitle = paste0("By SDI quartile. RMSE (", sample,")= ", formatC(rmse(preds_loc$population, preds_loc$pred_pop_cascade), digits = 4),
                           ". SDI coefficient = ", formatC(get(standard)$beta_soln[2], digits=3)),
         y = paste("Predicted ", measure),
         x = predictor,
         color = "SDI") +
    theme(legend.position="bottom",
          plot.title = element_text(size=8),
          plot.subtitle = element_text(size=7),
          legend.title.align = .5,
          legend.background = element_rect(fill="grey90", size=0.5, linetype="solid", colour ="grey30"),
          legend.key = element_rect(fill = "grey90")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    geom_vline(xintercept = c(get(standard)$cov_models[[3]]$spline_knots) ) +
    coord_cartesian(ylim=c(0, max(preds_loc$population, preds_loc$pred_pop, preds_loc$pred_pop_cascade)))
  
  
  
  gg_obs_spline <- ggplot(preds_spline, aes(y=pred_pop, x=get(exposure)))+
    geom_line(data = preds_spline,
              aes(color= sdi, group = sdi_group), alpha = 1, size=1)+
    theme_bw() +
    geom_line(data = preds_loc, aes(y=population, x=get(exposure), color= sdi, group = location_id), alpha = 0.7)+
    labs(title = title_obs_spline,
         subtitle = paste0("By SDI quartile. Observed vs. ", predictor),
         y = paste("Observed ", measure),
         x = predictor,
         color = "SDI") +
    theme(legend.position="bottom",
          plot.title = element_text(size=8),
          plot.subtitle = element_text(size=7),
          legend.title.align = .5,
          legend.background = element_rect(fill="grey90", size=0.5, linetype="solid", colour ="grey30"),
          legend.key = element_rect(fill = "grey90")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    geom_vline(xintercept = c(get(standard)$cov_models[[3]]$spline_knots) ) +
    coord_cartesian(ylim=c(0, max(preds_loc$population, preds_loc$pred_pop, preds_loc$pred_pop_cascade)))
  
  
  gg_stand_facet <- ggplot(preds_loc, aes(y=population, x=pred_pop))+
    geom_line(data = preds_loc,
              aes(color= region_name, group = location_id))+
    theme_bw() +
    scale_color_manual(values = custom.col.r)+
    labs(title = title_label_stand,
         subtitle = paste0("RMSE (", sample,")= ", formatC(rmse(preds_loc$population, preds_loc$pred_pop), digits = 4),
                           ". SDI coefficient = ", formatC(get(standard)$beta_soln[2], digits=3)),
         y = paste("Observed ", measure),
         x = paste("Predicted ", measure),
         color = "Super Region") +
    facet_wrap(~super_region_name)+
    theme(legend.position="none",
          plot.title = element_text(size=8),
          plot.subtitle = element_text(size=7),
          legend.title.align = .5,
          legend.background = element_rect(fill="grey90", size=0.5, linetype="solid", colour ="grey30"),
          legend.key = element_rect(fill = "grey90")) +
    guides(color = guide_legend(override.aes = list(alpha = .5), title.position = "top")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(ylim=c(0, max(preds_loc$population, preds_loc$pred_pop, preds_loc$pred_pop_cascade)), xlim=c(0, max(preds_loc$population, preds_loc$pred_pop, preds_loc$pred_pop_cascade)))+
    #identify the largest 5 residuals from the final year predicted
    geom_text(data=subset(preds_loc, abs(residual)%in%tail(sort( preds_loc[, max(abs(residual)), by = "location_id"]$V1 ),5)) , aes(y=population,x=pred_pop,label=paste0(location_name, " ", formatC(residual, digits=3, format="f") )), nudge_y = space,  size=1.3)
  
  gg_casc <- ggplot(preds_loc, aes(y=population, x=pred_pop_cascade))+
    geom_line(data = preds_loc,
              aes(color= super_region_name, group = location_id), alpha = 0.3)+
    theme_bw() +
    scale_color_manual(values = custom.col.sr)+
    labs(title = title_label_casc,
         subtitle = paste0("RMSE (", sample,")= ", formatC(rmse(preds_loc$population, preds_loc$pred_pop_cascade), digits = 4)),
         y = paste("Observed ", measure),
         x = paste("Predicted ", measure),
         color = "Super Region") +
    theme(legend.position="bottom",
          plot.title = element_text(size=8),
          plot.subtitle = element_text(size=7),
          legend.title.align = .5,
          legend.background = element_rect(fill="grey90", size=0.5, linetype="solid", colour ="grey30"),
          legend.key = element_rect(fill = "grey90")) +
    guides(color = guide_legend(override.aes = list(alpha = .5), title.position = "top")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(ylim=c(0, max(preds_loc$population, preds_loc$pred_pop, preds_loc$pred_pop_cascade)), xlim=c(0, max(preds_loc$population, preds_loc$pred_pop, preds_loc$pred_pop_cascade)))+
    #identify the largest 5 residuals from the final year predicted
    geom_text(data=subset(preds_loc, abs(residual_cascade) %in% tail(sort( preds_loc[, max(abs(residual_cascade)), by = "location_id"]$V1 ),5))  , aes(y=population,x=pred_pop_cascade,label=paste0(location_name, " ", formatC(residual_cascade, digits=3, format="f") )), nudge_y = space,  size=1.5)
  
  
  gg_casc_facet <- ggplot(preds_loc, aes(y=population, x=pred_pop_cascade))+
    geom_line(data = preds_loc,
              aes(color= region_name, group = location_id))+
    theme_bw() +
    scale_color_manual(values = custom.col.r)+
    labs(title = title_label_casc,
         subtitle = paste0("RMSE (", sample,")= ", formatC(rmse(preds_loc$population, preds_loc$pred_pop_cascade), digits = 4)),
         y = paste("Observed ", measure),
         x = paste("Predicted ", measure),
         color = "Super Region") +
    facet_wrap(~super_region_name)+
    theme(legend.position="none",
          plot.title = element_text(size=8),
          plot.subtitle = element_text(size=7),
          legend.title.align = .5,
          legend.background = element_rect(fill="grey90", size=0.5, linetype="solid", colour ="grey30"),
          legend.key = element_rect(fill = "grey90")) +
    guides(color = guide_legend(override.aes = list(alpha = .5), title.position = "top")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(ylim=c(0, max(preds_loc$population, preds_loc$pred_pop, preds_loc$pred_pop_cascade)), xlim=c(0, max(preds_loc$population, preds_loc$pred_pop, preds_loc$pred_pop_cascade)))+
    #identify the largest 5 residuals from the final year predicted
    geom_text(data=subset(preds_loc, abs(residual_cascade) %in% tail(sort( preds_loc[, max(abs(residual_cascade)), by = "location_id"]$V1 ),5))  , aes(y=population,x=pred_pop_cascade,label=paste0(location_name, " ", formatC(residual_cascade, digits=3, format="f") )), nudge_y = space,  size=1.3)
  #diagnostic type plots...time vs. SEV/population
  res_stand_facet <- ggplot(preds_loc, aes(y=residual, x=year_id))+
    geom_line(data = preds_loc,
              aes(color= region_name, group = location_id))+
    theme_bw() +
    scale_color_manual(values = custom.col.r)+
    labs(title = title_label_res_stand,
         subtitle = paste0("RMSE (", sample,")= ", formatC(rmse(preds_loc$population, preds_loc$pred_pop), digits = 4)),
         y = paste("Residuals"),
         x = paste("Year"),
         color = "Region") +
    facet_wrap(~super_region_name)+
    theme(legend.position="none",
          plot.title = element_text(size=8),
          plot.subtitle = element_text(size=7),
          legend.title.align = .5,
          legend.background = element_rect(fill="grey90", size=0.5, linetype="solid", colour ="grey30"),
          legend.key = element_rect(fill = "grey90")) +
    guides(color = guide_legend(override.aes = list(alpha = .5), title.position = "top")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    geom_abline(intercept = 0, slope = 0) +
    #identify the largest 5 residuals (by location) from across period
    geom_text(data=subset(preds_loc, abs(residual) %in% tail(sort( preds_loc[, max(abs(residual)), by = "location_id"]$V1 ),5)) , aes(y=residual,x=year_id,label=paste0(location_name, " ", formatC(residual, digits=3, format="f") )), hjust = 0, nudge_x = space_res , nudge_y = 0.001,  size=1.3)
  
  res_casc_facet <- ggplot(preds_loc, aes(y=residual_cascade, x=year_id))+
    geom_line(data = preds_loc,
              aes(color= region_name, group = location_id))+
    theme_bw() +
    scale_color_manual(values = custom.col.r)+
    labs(title = title_label_res_casc,
         subtitle = paste0("RMSE (", sample,")= ", formatC(rmse(preds_loc$population, preds_loc$pred_pop_cascade), digits = 4)),
         y = paste("Residuals"),
         x = paste("Year"),
         color = "Region") +
    facet_wrap(~super_region_name)+
    theme(legend.position="none",
          plot.title = element_text(size=8),
          plot.subtitle = element_text(size=7),
          legend.title.align = .5,
          legend.background = element_rect(fill="grey90", size=0.5, linetype="solid", colour ="grey30"),
          legend.key = element_rect(fill = "grey90")) +
    guides(color = guide_legend(override.aes = list(alpha = .5), title.position = "top")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    geom_abline(intercept = 0, slope = 0) +
    #identify the largest 5 residuals from the final year predicted
    geom_text(data=subset(preds_loc, abs(residual_cascade) %in% tail(sort( preds_loc[, max(abs(residual_cascade)), by = "location_id"]$V1 ),5))  , aes(y=residual_cascade,x=year_id,label=paste0(location_name, " ", formatC(residual_cascade, digits=3, format="f") )), hjust = 0, nudge_x = space_res , nudge_y = 0.001,  size=1.3)

  combined <- gg_obs_spline + gg_casc_spline_knots & theme(legend.position = "bottom")
  
  print(gg_stand)
  print(gg_stand_facet)
  print(res_stand_facet)
  print(gg_casc)
  print(gg_casc_facet)
  print(res_casc_facet)
}
}
dev.off()



#forecasts 2021-2050

cl <- makeCluster(getOption("cl.cores", 4))

forecasts_sex_population388<-do.call("rbind",parLapply(cl, c(1,2), function(sid){
  
  library(data.table)
  #MR-BRT and cascading splines functions require docker. Instructions are available from: https://github.com/ihmeuw-msca/mrtoolr/tree/main
  library(reticulate); use_condaenv("mrtool-0.0.1"); library(mrtoolr)
  mr <- import("mrtoolr")
  library(boot)
  library(zoo)
  invisible(sapply(list.files("FILEPATH", full.names = T), source))
  loc.met <- get_location_metadata(location_set_id = 92, release_id = 9)
  #all countries and territories
  loc.met <- loc.met[level==3]
  locs <- unique(loc.met$location_id)
  hierarchy <- get_location_metadata(92, release_id = 9)
  
  output_dir <- file.path("FILEPATH")
  output_dir2<-"FILEPATH"
  agid<- 388
  sdi_forecast <- fread("sdi_forecast.csv")# SDI forecasts
  sdi_forecast<- sdi_forecast[location_id%in%locs & year_id%in%c(2021:2050) & scenario== 0,]
  setnames(sdi_forecast, "mean", "sdi") 
  pop_forecast<- fread("population_forecast.csv")
  setnames(pop_forecast, old="value", new= "population")
  pop_2021 <- fread("population_retrospective.csv")
  pop_2021 <- pop_2021[year_id==2021 & age_group_id%in%c(1,2,3,4,5,8:14) & location_id!=1,]
  pop_forecast<- rbind(pop_forecast, pop_2021[,-"run_id"], fill= TRUE)
  pop_forecast<- pop_forecast[order(location_id, age_group_id, sex_id, year_id)]

  pop_forecast[age_group_id==3, pop3:= population]
  pop_forecast[, pop3:= na.locf(pop3, na.rm=FALSE), by= c("location_id", "sex_id", "year_id")]
  
  pop_forecast3_4 <- pop_forecast[location_id %in% locs & age_group_id==4]
  setnames(pop_forecast3_4, old="population", new= "pop4")
  pop_forecast3_4 <- merge(pop_forecast3_4, hierarchy[, .(location_id, super_region_name, region_name)], all.x=T, by = "location_id")
  pop_forecast3_4<- merge(pop_forecast3_4, sdi_forecast[scenario==0], by=c("location_id", "year_id"))
  
  pop_forecast3_4[,age_group_id:=388]
  pop_forecast3_4[,age_group_name := "1-5 Months"]
  standard <- "population_full"
  cascade <- "population_cascade_full"
  pop_388_df <- readRDS(file.path(output_dir, "ag388_population2021.RDS"))
  
  fit_pop_full <- function(model_data,sid) {
    train_df<- model_data[year_id>2004 & sex_id==sid,]

    #training data
    dat_loc <- MRData()
    dat_loc$load_df(
      data = train_df,
      col_obs = "population", col_obs_se = "se",
      col_covs = list("sdi", "pop3","pop4"), col_study_id = "location_id"
    )
    
    #"standard MRBRT model
    stand_mod <- MRBRT(
      data = dat_loc,
      cov_models = list(
        LinearCovModel("intercept", use_re = TRUE),
        LinearCovModel("sdi", use_re = FALSE),
        LinearCovModel("pop3", use_re = TRUE),
        LinearCovModel("pop4", use_re = FALSE, 
                       use_spline = TRUE,
                       spline_knots = array(seq(0, 1, by = 0.25)), # this tells is to put a spline every 0.25
                       spline_degree = 1L,
                       spline_knots_type = 'frequency',# this specifies to put the knots literally along the domain. Other option is frequency for it to be data density driven
                       spline_r_linear = TRUE,
                       spline_l_linear = FALSE
        )
      ))
    stand_mod$fit_model(inner_print_level = 0L, inner_max_iter = 1000L,  outer_max_iter= 500L)
    estimated_beta <- stand_mod$summary()[[1]][1, "sdi"]
    stand_mod$cov_models[[which(stand_mod$cov_names == "sdi")]]$prior_beta_uniform <- matrix(rep(estimated_beta, 2), ncol = 1)
    
    return(stand_mod)
    
  }
  
  fit_pop_cascade_full <- function(mod_global, model_data, output_dir, sid) {
    
    train_df<- model_data[sex_id==sid,]

    model_label_tmp <- paste0("cascade_population", unique(model_data$age_group_id), "_", sid)
    
    thetas <- c(2,7)
    cascade_fit <- run_spline_cascade(
      stage1_model_object = mod_global, 
      df = train_df, 
      col_obs = "population",
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
  
  
  population_full <- fit_pop_full(pop_388_df,sid)
  population_cascade_full <- fit_pop_cascade_full(population_full,pop_388_df,output_dir2,sid)
  
  forecast_df <- pop_forecast3_4[sex_id== sid]
  
  #prediction data- forecast
  dat_pred_forecast <- MRData()
  dat_pred_forecast$load_df(
    data = forecast_df,
    col_covs=list("sdi", "pop3","pop4"), col_study_id = "location_id")
  
  #Forecasting the future
  pred_pop <- get(standard)$predict(data = dat_pred_forecast, predict_for_study= TRUE, sort_by_data_id= TRUE)
  preds_loc_forecast <- forecast_df
  preds_loc_forecast$pred_pop <- pred_pop
  
  #forecasting the future
  preds_loc_forecast <- as.data.table(predict_spline_cascade(fit = get(cascade), newdata = preds_loc_forecast))
  preds_loc_forecast[, pred_pop_cascade := pred]

  return(preds_loc_forecast)
}))
forecasts_sex_population388[,location_id:= as.numeric(location_id)]
forecasts_sex_population388 <- merge(forecasts_sex_population388, hierarchy[, .(location_id, location_name)], all.x=T, by = "location_id")
dt<- forecasts_sex_population388[, c("location_id",  "age_group_id", "sex_id", "year_id", "pred_pop_cascade")]
setnames(dt, old= "pred_pop_cascade", new= "value")

#create sex_id 3 for pop_forecast 388
agg_dt<- dt
#sum the populations from the two sexes
agg_dt[, tot_pop_sex:= sum(value), by = c("location_id", "year_id", "age_group_id")]
agg_dt <- agg_dt[sex_id==1,-"value"]
agg_dt[,sex_id:=3]

setnames(agg_dt, old= "tot_pop_sex", new= "value")
#rbind back in the combined sex group
dt<- rbind(dt[,-"tot_pop_sex"], agg_dt)

saveRDS(dt, paste0("FILEPATH", "/pop388_locs", gsub("-", "", Sys.Date()) , ".RDS")) 
  