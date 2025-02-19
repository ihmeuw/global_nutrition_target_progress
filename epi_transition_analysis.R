################################################################################
## DESCRIPTION: Epidemiological Transition Analysis. Use meta-regression Bayesian prior tool (MR-BRT) to model the relationship between sociodemographic index (SDI) and prevalence of each GNT indicators. 
## Then use model to predict expected prevalence of each indicator from 2012 to 2021.
## INPUTS: 
#a) Separate country-age-sex-year prevalence estimates from 1990 to 2021 for Low-birthweight (LBW), exclusive breastfeeding (EBF), stunting, wasting, overweight, and anemia.
#b) Retrospective estimates of country-level socio-demographic index from 1990 to 2021.
## OUTPUTS: 
#a) Predicted prevalence for each GNT indicator from 2012 to 2021, root-mean squared error (RMSE)
################################################################################
#------------------SET-UP--------------------------------------------------
invisible(sapply(list.files("FILEPATH", full.names = T), source)) #IHME internal functions
path <- paste0("FILEPATH")

package_lib <- sprintf(path)
.libPaths(path)

#Color scheme for plots
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
# load packages, install if missing
library(reticulate)
library(scales)
packages <- c("data.table","magrittr","ggplot2", "Metrics", "boot", "tidyselect", "reticulate", "zoo")
for(p in packages){
  if(p %in% rownames(installed.packages())==FALSE){
    install.packages(p)
  }
  library(p, character.only = T)
}

#MR-BRT and cascading splines functions require docker. Instructions are available from: https://github.com/ihmeuw-msca/mrtoolr/tree/main
library(reticulate); use_condaenv("mrtool-0.0.1"); library(mrtoolr)
mr <- import("mrtoolr")

library(ModelMetrics)
library(data.table)
library(boot)
library(zoo)
# Directories -------------------------------------------------------------
out_dir <- file.path("FILEPATH")
invisible(sapply(list.files("/share/cc_resources/libraries/current/r/", full.names = T), source))
age_weights <- get_age_weights(release_id = 9)
hierarchy <- get_location_metadata(35, release_id = 9)
logit <- function(x){log(x/(1-x))}
invlogit <- function(x){exp(x)/(1+exp(x))}
prepare_df_format <- function(df, standard.locs){
  
  df[age_group_id == 1, age_group_name := "Under 5"]
  df[age_group_id == 2, age_group_name := "Early Neonatal"]
  df[age_group_id == 3, age_group_name := "Late Neonatal"]
  #for GBD2019 we need to use agegroupid 5 for overweight in chidren
  df[age_group_id == 5, age_group_name := "1 to 4"]
  df[age_group_id == 388, age_group_name := "1-5 Months"]
  df[age_group_id == 390, age_group_name := "<6 months"]
  df[age_group_id == 389, age_group_name := "6-11 Months"]
  df[age_group_id == 238, age_group_name := "12-23 Months"]
  df[age_group_id == 34, age_group_name := "2-4 Years"]
  df[age_group_id == 24, age_group_name := "15-49 years"]
  df[age_group_id == 164, age_group_name := "Birth"]
  
  df[modelable_entity_id==10556, indicator:= "stunting"]
  df[modelable_entity_id==10558, indicator:= "wasting"]
  df[modelable_entity_id==20018, indicator:= "overweight"]
  df[modelable_entity_id==16282, indicator:= "LBW"]
  df[modelable_entity_id==10507, indicator:= "anemia"]
  df[modelable_entity_id==20417, indicator:= "EBF"]
  
  loc.meta <- standard.locs[, c("location_id", "location_name", "location_type")]
  
  df <- merge(df, loc.meta, by = "location_id")
  
  df[measure_id == 5, measure_name := "Prevalence"]
  df[metric_id == 3, metric_name := "Rate"]
  df[sex_id == 1, sex := "Male"]
  df[sex_id == 2, sex := "Female"]
  df[, lyas := paste0(location_id, "_", year_id, "_", age_group_id, "_", sex_id)]
  df[, val := mean(value), by = lyas]
  df[, lower:= quantile(value, probs = .025), by = lyas]
  df[, upper:= quantile(value, probs = .975), by = lyas]
  df <- df[variable == "draw_0"]
  df$variable <- NULL
  df$value <- NULL
  
  return(df)
  
  
  
}
# Read in stunting draws and calculate means --------------------------------
id.vars <- c("metric_id", "age_group_id", "location_id", "measure_id", "modelable_entity_id", "sex_id", "year_id", "model_version_id")
standard.locs <- get_location_metadata(release_id = 9, location_set_id = 101)
standard.locations <- unique(standard.locs$location_id)
loc.met <- get_location_metadata(location_set_id = 91, release_id = 9)

#load formatted prvalence dataframe
epi_trans_df<-fread("FILEPATH", "retrospective_prevalence.csv")
epi_trans_df[modelable_entity_id==10556, indicator:= "stunting"]
epi_trans_df[modelable_entity_id==10558, indicator:= "wasting"]
epi_trans_df[modelable_entity_id==20018, indicator:= "overweight"]
epi_trans_df[modelable_entity_id==16282, indicator:= "LBW"]
epi_trans_df[modelable_entity_id==10507, indicator:= "anemia"]
epi_trans_df[modelable_entity_id==20417, indicator:= "EBF"]

epi_trans_df<- epi_trans_df[,metric_name:= "Rate"]
epi_trans_df<- epi_trans_df[indicator== "EBF",metric_name:= ""]

epi_trans_df<- epi_trans_df[indicator== "LBW",measure_name:= "Incidence"]
epi_trans_df<- epi_trans_df[indicator== "EBF",measure_name:= "Proportion"]
#proportion isn't a rate
epi_trans_df<- epi_trans_df[!(indicator== "LBW" | indicator== "EBF"),measure_name:= "Prevalence"]

# merge on sdi for locations
sdi <- get_covariate_estimates(881, release_id = 9)
epi_trans_df <- merge(epi_trans_df, sdi[,.(location_id, year_id, sdi = mean_value)], all.x=T, by = c("location_id", "year_id"))
epi_trans_df <- merge(epi_trans_df, hierarchy[, .(location_id, super_region_name, region_name)], all.x=T, by = "location_id")
epi_trans_df <- epi_trans_df[sex_id!=3,]
epi_trans_df$super_region_name <- factor(epi_trans_df$super_region_name)
epi_trans_df$region_name <- factor(epi_trans_df$region_name)

#calculate spearman correlation coefficients + significance for paper:
#######
cor_table<- data.frame(nrow=11, ncol=5)
i<-1
for(ind in c("LBW", "EBF", "stunting","wasting", "overweight", "anemia")){
  age <- min(epi_trans_df[indicator== ind,]$age_group_id)
  sid_list<- c("Male", "Female")
  if(ind== "anemia"){sid_list<- "Female"}
  for(sid in sid_list){
  test<- cor.test(epi_trans_df[indicator== ind & age_group_id== age & sex== sid,]$val, epi_trans_df[indicator==ind & age_group_id== age & sex== sid,]$sdi, method = "spearman")
  cor_table[i,1]<- ind
  cor_table[i,2]<- sid
  cor_table[i,3]<- formatC(test$estimate, digits = 4)
  cor_table[i, 4]<- formatC(test$p.value, digits = 4)
  cor_table[i, 5]<- test$statistic
  i<- i+1
  }

}
colnames(cor_table)<- c("Indicator", "Sex", "Spearman's rank rho", "p-value", "S" )
cor_table<- data.table(cor_table)
cor_table_g<- gt(cor_table[,c("Indicator", "Sex", "Spearman's rank rho", "p-value")])
write.csv(cor_table[,c("Indicator", "Sex", "Spearman's rank rho", "p-value")], paste0("FILEPATH"))


pdf(file=paste0("FILEPATH"), height = 11, width = 8.5)
cor_table_g<- tableGrob(cor_table[,c("Indicator", "Sex", "Spearman's rank rho")],
                               theme = ttheme_default(base_size = 10, base_family = "Times"), rows=c())

title <- textGrob("  Table S2: Spearman's rank correlation coefficients,\n by sex, standard GBD locations, 1990 - 2021", gp=gpar(fontsize=10, fontfamily = "Times"))
padding <- unit(2,"mm")

table <- gtable_add_rows(
  cor_table_g, 
  heights = grobHeight(title) + padding,
  pos = 0)
table <- gtable_add_grob(
  table, 
  title, 
  1, 1, 1, ncol(table))


grid.newpage()
grid.draw(table)
dev.off()

###################
# run age-, sex-, measure-specific splines
pdf(file.path(out_dir, paste0("GNT_spline_plots_", gsub("-", "", Sys.Date()),".pdf")), height = 7.5, width = 10)
out_dt <- data.table()
#for loop to run and store ggplots
for(ind in c("stunting", "anemia", "LBW", "EBF", "overweight", "wasting")){
  
  for(sid in unique(epi_trans_df[indicator==ind]$sex_id)){
    group <- "location_id"
      for(agid in unique(epi_trans_df[indicator==ind]$age_group_id)){
      data_dt <- epi_trans_df[indicator == ind & sex_id == sid & age_group_id == agid]
      # replace SE with the same value for all locations
      data_dt[, se := 1]
      
      # use offset before logit
      data_dt[val<1e-7, val := 1e-7]
      data_dt[val>(1-1e-7), val := (1-1e-7)]
      data_dt[, logit_val := logit(val)]
      

      # run MRBRT spline
      dat_1 <- MRData()
      dat_1$load_df(
        data = data_dt,  col_obs = "logit_val", col_obs_se = "se",
        col_covs = list("sdi"), col_study_id = group
      )
      
      if(ind== "overweight"){
        model <- MRBRT(data = dat_1,
                       cov_models = list(LinearCovModel("sdi", 
                                                        use_re = FALSE,
                                                        use_spline = TRUE,
                                                        spline_degree = 3L,
                                                        spline_knots = array(seq(0,1, by = .2)),
                                                        spline_knots_type = "frequency",
                                                        prior_spline_monotonicity= "increasing",
                                                        prior_spline_convexity = "convex",
                                                        spline_r_linear = FALSE,
                                                        spline_l_linear = TRUE),
                                         LinearCovModel("intercept", use_re = FALSE)), inlier_pct = .95)
      }
      
      if(ind == "EBF"){
        model <- MRBRT(data = dat_1,
                       cov_models = list(LinearCovModel("sdi", 
                                                        use_re = FALSE,
                                                        use_spline = TRUE,
                                                        spline_degree = 3L,
                                                        spline_knots = array(seq(0,1, by = .2)),
                                                        spline_knots_type = "frequency",
                                                        prior_spline_monotonicity= "increasing",
                                                        prior_spline_convexity = "concave",
                                                        spline_r_linear = FALSE,
                                                        spline_l_linear = TRUE),
                                         LinearCovModel("intercept", use_re = FALSE)), inlier_pct = .95)
      }
      
      if(!(ind== "overweight" | ind == "EBF")){
        model <- MRBRT(data = dat_1,
                       cov_models = list(LinearCovModel("sdi", 
                                                        use_re = FALSE,
                                                        use_spline = TRUE,
                                                        spline_degree = 3L,
                                                        spline_knots = array(seq(0,1, by = .25)),
                                                        spline_knots_type = "frequency",
                                                        prior_spline_monotonicity= "decreasing",
                                                        spline_r_linear = FALSE,
                                                        spline_l_linear = TRUE),
                                            LinearCovModel("intercept", use_re = FALSE)), inlier_pct = .95)
      }
      model$fit_model(inner_print_level = 5L, inner_max_iter = 200L, outer_max_iter = 500L)
      
      pred_dt <- rbind(data_dt, data.table(sdi = c(0,1)), fill = T)
      
      pred_dt <- pred_dt[order(sdi)]
      
      dat_pred <- MRData()
      
      dat_pred$load_df(
        data = pred_dt,
        col_covs=list("sdi")
      )
      
      pred_dt$expected <- model$predict(dat_pred) %>% invlogit

      out_dt <- rbind(out_dt, pred_dt[!is.na(location_id)])
      
      print(paste("Done with", ind, sid, agid, group))
    }
    
    out_dt$age_group_name <- factor(out_dt$age_group_name, levels = c("Birth","Early Neonatal", "Late Neonatal", "1-5 Months", "2-4 Years", "Under 5", "15-49 years"))
    out_dt$super_region_name <- factor(out_dt$super_region_name)
    out_dt$region_name <- factor(out_dt$region_name)
    
    for(agid in unique(epi_trans_df[indicator==ind]$age_group_id)){
      #main plot with all level 3 locations and global spline
      gg <- ggplot(out_dt[sex_id == sid & indicator==ind & age_group_id == agid], aes(y=val, x=sdi))+
        geom_line(data = epi_trans_df[sex_id == sid & indicator==ind  & age_group_id == agid & year_id %in% c(2012:2021)], 
                  aes(color = super_region_name, group = location_id), alpha = 0.2)+
        geom_line(aes(y = expected))+
        theme_bw()+
        scale_color_manual(values = custom.col.sr)+
        labs(title = paste0("Overall ", ind, " ", unique(data_dt$measure_name), " among ", 
                            unique(data_dt$sex), "s, ", unique(out_dt[indicator==ind & age_group_id == agid]$age_group_name),", 2012-2021, all locations"),
             subtitle = paste0("RMSE (2012-2021)= ", formatC(rmse(out_dt[sex_id == sid & indicator==ind & age_group_id == agid & year_id %in% c(2012:2021)]$val, out_dt[sex_id == sid & indicator==ind & age_group_id == agid & year_id %in% c(2012:2021)]$expected), digits = 3), "."),# groupvar= ", group),
              y = paste(unique(data_dt$measure_name), unique(data_dt$metric_name)), 
             x = "Sociodemographic Index",
             color = "Super Region")+
        theme(legend.position= "bottom", 
              legend.title.align = .5,
              legend.background = element_rect(fill= "grey90", size=0.5, linetype= "solid", colour = "grey30"),
              legend.key = element_rect(fill = "grey90")) +
        guides(color = guide_legend(override.aes = list(alpha = .5), title.position = "top")) +
        scale_y_continuous(labels = percent_format(accuracy = 1)) +
        theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
      
      
      print(gg)
      
    }
    }
  }


dev.off()

#what is the RMSE from the epi transition analyses?
for(ind in c("stunting", "anemia", "LBW", "EBF", "overweight", "wasting")){
  for(sid in unique(out_dt[indicator==ind]$sex_id)){
    for(agid in unique(out_dt[indicator==ind]$age_group_id)){
      print(paste0(ind, " ", unique(out_dt[indicator==ind & sex_id==sid]$sex)," , " ,agid, ", RMSE= ", rmse(out_dt[indicator==ind & sex_id==sid & age_group_id==agid]$val, out_dt[indicator==ind & sex_id==sid & age_group_id==agid]$expected)))
    }
  }
}


#save without the "both" sex_id category
write.csv(out_dt[sex_id!=3,], file.path(out_dir, "observed_and_expected_se1_v4.csv"), row.names = F)
out_dt <- fread(file.path(out_dir, "observed_and_expected_se1_v4.csv"))
out_dt <- out_dt[, .(indicator, location_id, age_group_id, age_group_name, year_id, sex_id, sex, measure_id, metric_id, sdi, expected, observed = val, measure_name, metric_name)]

# Create aggregates -------------------------------------------------------
# merge on populations
population <- get_population(release_id = 9,
                             location_id = unique(out_dt$location_id), 
                             age_group_id = unique(out_dt$age_group_id),
                             sex_id = unique(out_dt$sex_id),
                             year_id = unique(out_dt$year_id))
population[, run_id := NULL]
out_dt <- merge(out_dt, population, by = c("location_id", "age_group_id", "year_id", "sex_id"), all.x=T)
# merge on age weights
out_dt <- merge(out_dt, age_weights, by = "age_group_id", all.x=T)
# Aggregate up location hierarchy
for(loc in hierarchy[most_detailed == 0, location_id]){
  
  if(loc == 1){
    agg_dt <- copy(out_dt)
  }else{
    child_locs <- hierarchy[most_detailed == 1 & (parent_id == loc | grepl(paste0(",", loc, ","), path_to_top_parent)), location_id]
    
    agg_dt <- out_dt[location_id %in% child_locs] 
  }
  
  agg_dt <- agg_dt[, .(location_id = loc, expected = weighted.mean(expected, w = population), observed = weighted.mean(observed, w = population),
                       population = sum(population)),
                   by = c("indicator","age_group_id", "age_group_name", "year_id", "sex_id", "sex","measure_id", "metric_id", "measure_name", "metric_name", "sdi","age_group_weight_value")] #"group_var"
  
  out_dt <- rbind(out_dt, agg_dt, use.names = T)
  
}
# Aggregate to desired age groups
ages <- list(
            "Birth" = 164,
            "<6 months" = c(2, 3, 388),
             "2-4 years" = 34,
             "<5 years" = 1,
             "15-49 years" = 24,
             "All Ages" = unique(out_dt$age_group_id),
             "Age-standardized" = unique(out_dt$age_group_id))
sum_dt <- data.table()
for(i in 1:length(ages)){
  
  agg_dt <- out_dt[age_group_id %in% ages[[i]]]
  
  if(names(ages)[i]== "Age-standardized"){
    agg_dt[, weight := age_group_weight_value]
  }else{
    agg_dt[, weight := population]
  }
  
  agg_dt <- agg_dt[, .(expected = weighted.mean(expected, w = weight), observed = weighted.mean(observed, w = weight), age_group_name = names(ages)[i]),
                   by = c("location_id", "indicator","year_id", "measure_id", "metric_id", "measure_name",  "metric_name")]
  
  sum_dt <- rbind(sum_dt, agg_dt, use.names = T)
  
}
# merge on location names
sum_dt <- merge(sum_dt, hierarchy[, .(location_id, location_name, ihme_loc_id, level, super_region_name, region_name)], by = "location_id", all.x=T)
sum_dt[, ratio := observed/expected]
sum_dt <- merge(sum_dt, sdi[,.(location_id, year_id, sdi = mean_value)], all.x=T, by = c("location_id", "year_id"))

for(ind in c("stunting", "anemia", "LBW", "EBF", "overweight", "wasting")){
    for(agid in unique(sum_dt[indicator==ind & age_group_name!= "All Ages" & age_group_name!= "Age-standardized"]$age_group_name)){
      print(paste0(ind, ", ", agid,  ", RMSE= ", rmse(sum_dt[indicator==ind & level==3 & age_group_name==agid]$observed, sum_dt[indicator==ind & level==3 & age_group_name==agid]$expected)))
      
    }
  }

write.csv(sum_dt, file.path(out_dir, "summary_results_se1_2021.csv"), row.names = F)

sum_dt <- as.data.table(read.csv(paste0(out_dir, "summary_results_se1_2021.csv")))
sum_dt <- sum_dt[!(age_group_name== "Age-standardized" | age_group_name== "All Ages"),]
sum_dt[, measure_name:= "Prevalence"]
#multiply sdi by 100
sum_dt[,sdi:=sdi*100]
sum_dt[indicator== "stunting", label:= "Stunting"]
sum_dt[indicator== "wasting", label:= "Wasting"]
sum_dt[indicator== "LBW", label:= "Low birthweight"]
sum_dt[indicator== "EBF", label:= "Exclusive breastfeeding"]
sum_dt[indicator== "overweight", label:= "Overweight"]
sum_dt[indicator== "anemia", label:= "Anaemia"]

out_dt[,sdi:= sdi*100]

pdf(file.path(out_dir, paste0("GNT_spline_plots_agg_", gsub("-", "", Sys.Date()),".pdf")), height = 7.5, width = 10)

for(ind in c("stunting", "anemia", "LBW", "EBF", "overweight", "wasting")){
  for(agid in unique(sum_dt[indicator==ind & age_group_name!= "All Ages" & age_group_name!= "Age-standardized"]$age_group_name)){
    sex <- ifelse(ind== "anemia", "females", "both sexes")

    data_dt <- sum_dt[indicator == ind & age_group_name == agid]
    
    #getting 500 SDI values along this span to smooth out
    min.sdi <- min(out_dt$sdi)
    max.sdi <- max(out_dt$sdi)
    sdi.vals <- seq(from = min.sdi, to = max.sdi, length.out = 500) # A vector of length 500 evenly spaced from the lowest to highest sdi values
    subset <- out_dt[indicator==ind & age_group_id==max(out_dt[indicator==ind]$age_group_id),]
    #getting the closest spline values to those SDI points
    overall.sdi.spaced.vals <- lapply(sdi.vals, function(sval){ # Find the spline values (y values) that have corresponding x values closest to those 500 points
      
      estimate <- subset[which.min(abs(sval-sdi))]$expected # Store the y value for the dataset row which had the x value closest to the x value in the vector
      
      points <- data.table(sdi_value = sval, overall.spline = estimate)  #creates a data table with that x value from the vector, and the corresponding y value which was identified above
      
      return(points)
      
    }) %>% rbindlist() # We now have a 500 row dataframe with x and y values

    #Calculate a rolling average for every 10 values
    overall.sdi.spaced.vals <- overall.sdi.spaced.vals %>% 
      mutate(roll_mean = rollmean(overall.spline, 10, na.pad = T)) #create a new column called roll_mean which takes a rolling average of overall.spline, which was the y value
    
    overall.sdi.spaced.vals <- data.table(overall.sdi.spaced.vals)
    setnames(overall.sdi.spaced.vals, old= "sdi_value", new = "sdi")
    
    #main plot with all level 3 locations and global spline
  gg <- ggplot(sum_dt[level >= 3 & indicator==ind & age_group_name == agid], aes(y=observed, x=sdi))+
    geom_line(data = sum_dt[indicator==ind & level==3 & year_id %in% c(2012:2021)], 
              aes(color = super_region_name, group = location_id), alpha = 0.7, size=0.5)+
    geom_line(data = overall.sdi.spaced.vals, aes(y = roll_mean), size=0.75)+
    theme_bw()+
    scale_color_manual(values = custom.col.sr)+
    labs(title = paste0("Co-evolution of ", ind, " ", unique(data_dt$measure_name), " with SDI, ", 
                        sex, ", ", unique(sum_dt[indicator==ind & age_group_name == agid]$age_group_name),", 2012 to 2021, 204 countries and territories"),
         subtitle = paste0("RMSE (IS)= ", formatC(rmse(sum_dt[indicator==ind & level==3 & age_group_name==agid]$observed, sum_dt[indicator==ind & level==3 & age_group_name==agid]$expected), digits = 3), "."),
         y = paste(unique(data_dt$measure_name)), 
         x = "Sociodemographic Index",
         color = "Super Region")+
    theme(legend.position= "bottom", 
          legend.title.align = .5,
          legend.background = element_rect(fill= "grey90", size=0.5, linetype= "solid", colour = "grey30"),
          legend.key = element_rect(fill = "grey90")) +
    guides(color = guide_legend(override.aes = list(alpha = .5), title.position = "top")) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
  
  
  print(gg)
  
    }
  }

dev.off()


library(dplyr)
epi_trans_summary_results <- as.data.table(read.csv(file.path(out_dir, "summary_results_se1_2021.csv")))
epi_trans_age_group <- epi_trans_summary_results[!(age_group_name== "Age-standardized" | age_group_name== "All Ages"),]

#observed rate of change from 2012-2021 vs expected rate of change from 2012-2021
#AROC calculation since 2012
epi_trans_age_group[, AROC_expected:= log( expected / lag(expected,9) ) / 9, by = c("indicator", "location_id")]
epi_trans_age_group[year_id!=2021, AROC_expected:= NA ]

epi_trans_age_group[, AROC_observed:= log( observed / lag(observed,9) ) / 9, by = c("indicator", "location_id")]
epi_trans_age_group[year_id!=2021, AROC_observed:= NA ]
#relative difference in rates of change
epi_trans_age_group[year_id==2021, AROC_difference:= AROC_observed - AROC_expected]



#merge in the 80th percentile for indicator-specific ARC draws to give some more sense of uncertainty of the estimates
library(stringr)
for(ind in c("LBW", "EBF", "overweight", "stunting", "wasting", "anemia")){
  ind_name<- ind
  if(ind!= "LBW"){ind_name<- str_to_title(ind)}
  if(ind== "EBF"){ ind_name<- "Breastfeeding"}
  
  indicator_df<-readRDS(paste0("FILEPATH", ind_name, "/estimates1990_2021.RDS"))
  sid<-max(unique(indicator_df$sex_id))
  epi_trans_temp<-merge(epi_trans_age_group[indicator== ind & level %in% c(0:3) & year_id==2021,], indicator_df[sex_id == sid, c("location_id", "year_id", "arc_80pct") ], by =c("location_id", "year_id"), all.x= TRUE)
  #rbind the indicaor-specific data frames back together
  if(ind== "LBW"){
    epi_trans_all<- epi_trans_temp
  }else{
    epi_trans_all<- rbind(epi_trans_all, epi_trans_temp)
  }
  
}

#Save year 2021 epi-transition results
write.csv(epi_trans_all[year_id==2021 & level %in% c(0:3), -c("outperform_change")], file.path(out_dir, "epi_trans_compare_y2021.csv"), row.names = F)
#Save year 2012 epi-transition results
write.csv(epi_trans_age_group[year_id==2012 & level %in% c(0:3), -c("AROC_expected", "AROC_observed", "AROC_difference")], file.path(out_dir, "epi_trans_compare_y2012.csv"), row.names = F)

