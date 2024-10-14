################################################################################
## DESCRIPTION: Use meta-regression Bayesian prior tool (MR-BRT) and cascading splines to model prevalence from summary exposure values (SEVs)
## INPUTS: Low-birthweight (LBW), exclusive breastfeeding (EBF), stunting, wasting, overweight, and anemia SEVs, prevalence, and socio-demographic index (SDI) 
## OUTPUTS: Predicted prevalence in training sample (1990-2014), predicted prevalence in testing sample (2015-2021), root-mean squared error (RMSE) ##
## AUTHOR: 
## DATE:
################################################################################

#
library(tidyverse)
library(haven)
library(parallel)
invisible(sapply(list.files("FILEPATH", full.names = T), source)) #IHME internal functions

path <- paste0("FILEPATH")
#MR-BRT and cascading splines functions require docker. Instructions are available from: https://github.com/ihmeuw-msca/mrtoolr/tree/main
library(reticulate); use_condaenv("mrtool-0.0.1"); library(mrtoolr)
mr <- import("mrtoolr")

# Directories -------------------------------------------------------------
output_dir <- file.path("FILEPATH")

package_lib <- sprintf(path)
.libPaths(path)

library(matrixStats)
library(data.table)
library(zoo)
library(plyr)
library(dplyr)
library('reticulate')
library(scales)
library(ggplot2)
library(ModelMetrics)
library(boot)
packages <- c("magrittr","ggplot2", "Metrics", "boot", "tidyselect", "grid", "gridExtra", "ggpubr", "patchwork")
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

#Function for formatting data as desired
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
id.vars <- c("metric_id", "age_group_id", "location_id", "measure_id", "modelable_entity_id", "sex_id", "year_id", "model_version_id")
hierarchy <- get_location_metadata(92, release_id = 9)

#sdi
sdi <- fread("sdi_retrospective.csv")
setnames(sdi, old=c("mean","lower", "upper"), new=c("sdi","lower_sdi", "upper_sdi"))


# SEVs for all most detailed risk factors- GBD2021 versions
#CGF (wasting and stunting) create separate sevs for GBD2019 age_group_ids 4 (388, 389) and 5 (34, 238). These groups no longer exist in GBD 2021 but are used in forecasted SEVs.
sev_stunting = get_outputs("rei", rei_id=241, measure_id=29, metric_id=3, location_id = locs, age_group_id = 1, sex_id = c(1,2,3), release_id= 9, compare_version_id= 7982, year_id=c(1990:2021))
sev_stunting_45 = get_outputs("rei", rei_id=241, measure_id=29, metric_id=3, location_id = locs, age_group_id = c(238, 34, 388,389), sex_id = c(1,2,3), release_id= 9, compare_version_id= 7982 , year_id=c(1990:2021))

sev_wasting = get_outputs("rei", rei_id=240, measure_id=29, metric_id=3, location_id = locs, age_group_id = 1, sex_id = c(1,2,3), release_id= 9, compare_version_id= 7982, year_id=c(1990:2021))
sev_wasting_45 = get_outputs("rei", rei_id=240, measure_id=29, metric_id=3, location_id = locs, age_group_id = c(238, 34, 388,389), sex_id = c(1,2,3), release_id= 9, compare_version_id= 7982, year_id=c(1990:2021))

#simplified (combined agegroup sev for the diagnostics plots)
sev_irondef = get_outputs("rei", rei_id=95, measure_id=29, metric_id=3, location_id = locs, age_group_id = 24, sex_id = 2, release_id= 9, compare_version_id= 7999, year_id=c(1990:2021))
sev_ebf = get_outputs("rei", rei_id=136, measure_id=29, metric_id=3, location_id = locs, age_group_id = 3, sex_id = c(1,2,3), release_id= 9, compare_version_id= 7982, year_id=c(1990:2021))
sev_lbw = get_outputs("rei", rei_id=335, measure_id=29, metric_id=3, location_id = locs, age_group_id = 2, sex_id = c(1,2,3), release_id= 9, compare_version_id= 7982, year_id=c(1990:2021))
#overweight SEVs from children and adults
sev_overweight = get_outputs("rei", rei_id=371, rei_set_id=2, measure_id=29, metric_id=3, location_id = locs, age_group_id = 6, sex_id = c(1,2,3), release_id= 9, compare_version_id= 7982, year_id=c(1990:2021))

 #added adult overweight SEV back on August 26- updated compare_version_id to newest (Oct7) on Nov 14 then final on Jan 2 2024
 sev_overweight_adult = get_outputs("rei", rei_id=370, rei_set_id=2, measure_id=29, metric_id=3, location_id = locs, age_group_id = c(9:14), sex_id = c(1,2,3), release_id= 9, compare_version_id= 7982, year_id=c(1990:2021))
 
 populations <- get_population(age_group_id = c(9:14), location_set_id = 91, year_id = c(1990:2021), release_id = 9, location_id = locs, sex_id = c(1,2,3))
 sev_overweight_adult<- merge(sev_overweight_adult, populations, by= c("age_group_id", "location_id",  "year_id", "sex_id"))
 sev_overweight_adult[, sev:= weighted.mean(val, population), by= c("location_id", "sex_id", "year_id")]
 sev_overweight_adult<- sev_overweight_adult[age_group_id==9,]
 sev_overweight_adult[, val:= sev]
 
 setnames(sev_overweight, old=c("val","upper","lower"), new= c("sev_kids","upper_sev_kids","lower_sev_kids"))
 setnames(sev_overweight_adult, old=c("val","upper","lower"), new= c("sev_adults","upper_sev_adults","lower_sev_adults"))
 
sev_overweight_merge<-merge(sev_overweight, sev_overweight_adult[,c(2:4,18:20)], by= c("location_id", "year_id", "sex_id"))
sev_overweight_merge<- merge(sev_overweight_merge, sdi[,c("location_id", "year_id", "sdi", "lower_sdi", "upper_sdi")], by=c("location_id", "year_id"))

#stunting
#Prevalence estimates from GBD2021- include each of the neonatal agegroups
estimates1990_2021<-get_draws("modelable_entity_id", 10556, year_id=c(1990:2021),
                     source="epi", release_id = 9,
                     age_group_id = c(2,3,4,5), sex_id = c(1, 2, 3), location_id = locs)


estimates1990_2021 <- melt(estimates1990_2021, id.vars = id.vars)
estimates1990_2021 <- prepare_df_format(df = estimates1990_2021, loc.met[level==3])
stunt_prev <- merge(estimates1990_2021, hierarchy[, .(location_id, super_region_name, region_name)], all.x=T, by = "location_id")

stunt_prev <- stunt_prev[location_id%in%loc.met$location_id, c("location_id", "year_id", "age_group_id","sex_id","super_region_name","region_name","val", "lower","upper")]
setnames(stunt_prev, old=c("val","lower", "upper"), new=c("prevalence","lower_prev", "upper_prev"))

populations <- get_population(age_group_id = c(388,389), location_set_id = 91, year_id = c(1990:2021), release_id = 9, location_id = locs, sex_id = c(1,2,3))

#population-weight SEVs from agegroup 388 and 389 to get a SEV for agegroup 4 (post-neonatal)
sev_stunting_4<- sev_stunting_45[age_group_id%in%c(388,389),]
sev_stunting_4 <- merge(sev_stunting_4, populations[,-"run_id"], by= c("age_group_id", "location_id",  "year_id", "sex_id"))
sev_stunting_4[age_group_id>380, sev:= weighted.mean(val, population), by = c("location_id",  "measure_id",   "metric_id", 
                                                                              "rei_id", "sex_id", "year_id", "expected",
                                                                              "location_name", "location_type", "measure_name",  "metric_name",  
                                                                              "rei", "rei_name", "sex" )]
sev_stunting_4[age_group_id>380, lower_sev:= weighted.mean(lower, population), by = c("location_id",  "measure_id",   "metric_id", 
                                                                                      "rei_id", "sex_id", "year_id", "expected",
                                                                                      "location_name", "location_type", "measure_name",  "metric_name",  
                                                                                      "rei", "rei_name", "sex" )]
sev_stunting_4[age_group_id>380, upper_sev:= weighted.mean(upper, population), by = c("location_id",  "measure_id",   "metric_id", 
                                                                                      "rei_id", "sex_id", "year_id", "expected",
                                                                                      "location_name", "location_type", "measure_name",  "metric_name",  
                                                                                      "rei", "rei_name", "sex" )]
#drop val, lower, and upper
sev_stunting_4<- sev_stunting_4[age_group_id==388, -c("val","lower", "upper","population")]
sev_stunting_4[, age_group_id:=4]
sev_stunting_4[, age_group_name:= "Post Neonatal"]
#merge sev_stunting_4 with stunt_prev and then population-weight the prevalence to get prevalence for agegroup 4
stunt_df4<- merge(sev_stunting_4, stunt_prev, by=c("location_id", "year_id", "sex_id", "age_group_id"))

#fill the SEVs from age_group_id==4 into SEV for age_group_id==2 and 3
stunt_df23<-merge(sev_stunting_4[, -c("age_group_id","age_group_name")], stunt_prev[age_group_id%in%c(2,3)], by=c("location_id", "year_id", "sex_id"))

#population-weight SEVs from agegroup 238 and 34 to get a sev for agegroup 5 (1 to 4 years)
sev_stunting_5<- sev_stunting_45[age_group_id%in%c(238,34),]

populations <- get_population(age_group_id = c(238,34), location_set_id = 91, year_id = c(1990:2021), release_id = 9, location_id = locs, sex_id = c(1,2,3))

#population-weight SEVs from agegroup 388 and 389 to get a sev for agegroup 4 (post-neonatal)
sev_stunting_5 <- merge(sev_stunting_5, populations[,-"run_id"], by= c("age_group_id", "location_id",  "year_id", "sex_id"))
sev_stunting_5[, sev:= weighted.mean(val, population), by = c("location_id",  "measure_id",   "metric_id", 
                                                                              "rei_id", "sex_id", "year_id", "expected",
                                                                              "location_name", "location_type", "measure_name",  "metric_name",  
                                                                              "rei", "rei_name", "sex" )]
sev_stunting_5[, lower_sev:= weighted.mean(lower, population), by = c("location_id",  "measure_id",   "metric_id", 
                                                                                      "rei_id", "sex_id", "year_id", "expected",
                                                                                      "location_name", "location_type", "measure_name",  "metric_name",  
                                                                                      "rei", "rei_name", "sex" )]
sev_stunting_5[, upper_sev:= weighted.mean(upper, population), by = c("location_id",  "measure_id",   "metric_id", 
                                                                                      "rei_id", "sex_id", "year_id", "expected",
                                                                                      "location_name", "location_type", "measure_name",  "metric_name",  
                                                                                      "rei", "rei_name", "sex" )]
#select just one agegroup to sub in for age_group 5
sev_stunting_5<- sev_stunting_5[age_group_id==34, -c("val","lower", "upper", "population")]
sev_stunting_5[,age_group_id:=5]
sev_stunting_5[, age_group_name:= "1 to 4"]
#merge sev_stunting_4 with stunt_prev and then population-weight the prevalence to get prevalence for agegroup 4
stunt_df5<- merge(sev_stunting_5, stunt_prev, by=c("location_id", "year_id", "sex_id", "age_group_id"))

#rbind together stunt_df and stunt_df_23
stunt_df <- rbind(stunt_df23, stunt_df4, stunt_df5, fill=TRUE)
stunt_df[age_group_id == 2, age_group_name := "Early Neonatal"]
stunt_df[age_group_id == 3, age_group_name := "Late Neonatal"]

stunt_df<- merge(stunt_df, sdi[,c("location_id", "year_id", "sdi", "lower_sdi", "upper_sdi")], by=c("location_id", "year_id"))
#correlation for agegroup2 prevalence with agegroup4 SEV
cor(stunt_df[sex_id==2 & age_group_id==2]$prevalence, stunt_df[sex_id==2 & age_group_id==2]$sev)
cor(stunt_df[year_id==1990 & sex_id ==3]$prevalence, stunt_df[year_id==1990 & sex_id==3]$sev)
cor(stunt_df[year_id==2010 & sex_id ==3]$prevalence, stunt_df[year_id==2010 & sex_id==3]$sev)
cor(stunt_df[year_id==2019 & sex_id ==3]$prevalence, stunt_df[year_id==2019 & sex_id==3]$sev)

plot(stunt_df$sev, stunt_df$prevalence)
#plot for stunting
stunt_plot <- ggplot(stunt_df[sex_id==3 & age_group_id==5], aes(x=sev, y=prevalence, color=region_name)) +
  geom_point()+
  labs(title = "Stunting prevalence vs. SEV (1990- 2021), 1 to 4 years") +
  xlab("SEV") +
  ylab("Prevalence")+
  theme(legend.title = element_text(size = 4), legend.text = element_text(size = 4)) +
  geom_abline()+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  labs(subtitle=paste0("Correlation= ", formatC(cor(stunt_df[sex_id==3 & age_group_id==5]$prevalence, stunt_df[sex_id==3 & age_group_id==5]$sev),digits=3, format="f")) )


#irondef/anemia
 estimates1990_2021<-as.data.table(readRDS("FILEPATH"))
 estimates1990_2021<- estimates1990_2021[, -c("population", "run_id")]
 estimates1990_2021 <- melt(estimates1990_2021, id.vars = id.vars)
 estimates1990_2021 <- prepare_df_format(df = estimates1990_2021, loc.met[level==3])
 anem_prev <- merge(estimates1990_2021, hierarchy[, .(location_id, super_region_name, region_name)], all.x=T, by = "location_id")
 anem_prev <- anem_prev[location_id%in%loc.met$location_id & sex_id==2, c("location_id", "year_id", "age_group_id","super_region_name","region_name","val", "lower","upper")]
 setnames(anem_prev, old=c("val","lower", "upper"), new=c("prevalence","lower_prev", "upper_prev"))
 setnames(sev_irondef, old=c("val","lower", "upper"), new=c("sev","lower_sev", "upper_sev"), skip_absent=TRUE)

 anem_df <- merge(sev_irondef, anem_prev, by=c("location_id", "year_id", "age_group_id"))
cor(anem_df$prevalence, anem_df$sev)

anem_df<- merge(anem_df, sdi[,c("location_id", "year_id", "sdi", "lower_sdi", "upper_sdi")], by=c("location_id", "year_id"))


#Iron Deficiency SEV vs anemia prevalence:
anem_plot <- ggplot(anem_df, aes(x=sev, y=prevalence, color=region_name)) +
  geom_point()+
  labs(title = "Anemia prevalence vs. Iron deficiency SEV (1990-2021), females, 15-49 years") +
  xlab("SEV") +
  ylab("Prevalence")+
  theme(plot.title = element_text(size=9), 
        plot.subtitle = element_text(size=7), 
        legend.title = element_text(size = 4), 
        legend.text = element_text(size = 4),
        legend.title.align = .5,
        legend.background = element_rect(fill="grey90", size=0.5, linetype="solid", colour ="grey30"),
        legend.key = element_rect(fill = "grey90"),
        legend.position="right") +
  geom_abline()+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  labs(subtitle=paste0("Correlation= ", formatC(cor(anem_df$prevalence, anem_df$sev),digits=3, format="f")))


#LBW
LBW_prev <- readRDS("FILEPATH")
setnames(LBW_prev, old="mean", new="prevalence")

LBW_prev <- LBW_prev[location_id%in%loc.met$location_id, c("location_id", "year_id", "age_group_id","sex_id","super_region_name","region_name","prevalence", "lower_prev","upper_prev")]

setnames(sev_lbw, old=c("val","lower", "upper"), new=c("sev","lower_sev", "upper_sev"))

lbw_df<-merge(sev_lbw, LBW_prev[,-c("age_group_id")], by=c("location_id", "year_id","sex_id"))
lbw_df<- merge(lbw_df, sdi[,c("location_id", "year_id", "sdi", "lower_sdi", "upper_sdi")], by=c("location_id", "year_id"))
cor(lbw_df$prevalence, lbw_df$sev)
cor(lbw_df[year_id==1990]$prevalence, lbw_df[year_id==1990]$sev)
cor(lbw_df[year_id==2010]$prevalence, lbw_df[year_id==2010]$sev)
cor(lbw_df[year_id==2019]$prevalence, lbw_df[year_id==2019]$sev)

plot(lbw_df$sev, lbw_df$prevalence)

#plot for LBW
lbw_plot <- ggplot(lbw_df[year_id%in%c(1990, 1995, 2000, 2005, 2010, 2015, 2017, 2019)], aes(x=sev, y=prevalence, color=region_name)) +
  geom_point()+
  labs(title = "LBW incidence vs. SEV (1990, 1995, 2000, 2005, 2010, \n 2015, 2017, 2019)- GBD 2021") +
  xlab("SEV") +
  ylab("Incidence")+
  theme(legend.title = element_text(size = 4), legend.text = element_text(size = 4)) +
  geom_abline()+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  labs(subtitle=paste0("Correlation= ", formatC(cor(lbw_df[year_id%in%c(1990, 1995, 2000, 2005, 2010, 2015, 2017, 2019)]$prevalence, lbw_df[year_id%in%c(1990, 1995, 2000, 2005, 2010, 2015, 2017, 2019)]$sev),digits=3, format="f")) )


#Overweight
df<-readRDS("FILEPATH")
df<- df[location_id%in%loc.met[level==3]$location_id]
estimates1990_2021 <- melt(df, id.vars = id.vars)
estimates1990_2021 <- prepare_df_format(df = estimates1990_2021, loc.met[level==3])

overweight_prev <- merge(estimates1990_2021, hierarchy[, .(location_id, super_region_name, region_name)], all.x=T, by = "location_id")

overweight_prev <- overweight_prev[location_id%in%loc.met$location_id, c("location_id", "year_id", "sex_id","super_region_name","region_name","val", "lower","upper")]
setnames(overweight_prev, old=c("val","lower", "upper"), new=c("prevalence","lower_prev", "upper_prev"))
overweight_df<-merge(sev_overweight_merge[,-c("age_group_id")], overweight_prev, by=c("location_id", "year_id","sex_id"))
overweight_df[, age_group_name:= "2-4 Years"]

#EBF
estimates1990_2021<-as.data.table(readRDS("FILEPATH"))

#from GBD 2019, had to use late neonatal agegroup since old agegroup_set was used
estimates1990_2021<-get_draws("modelable_entity_id", 20417, year_id=c(1990:2021),
                              source="epi", release_id = 9,
                              age_group_id = c(2,3,388) , sex_id = c(1, 2, 3), location_id = locs)
estimates1990_2021 <- melt(estimates1990_2021, id.vars = id.vars)
estimates1990_2021 <- prepare_df_format(df = estimates1990_2021, loc.met[level==3])

ebf_prev <- merge(estimates1990_2021, hierarchy[, .(location_id, super_region_name, region_name)], all.x=T, by = "location_id")

ebf_prev <- ebf_prev[location_id%in%loc.met$location_id, c("location_id", "year_id", "sex_id", "age_group_id","super_region_name","region_name","val", "lower","upper")]
setnames(ebf_prev, old=c("val","lower", "upper"), new=c("prevalence","lower_prev", "upper_prev"))
setnames(sev_ebf, old=c("val","lower", "upper"), new=c("sev","lower_sev", "upper_sev"))

#We will use the same SEV to predict EBF prevalence in under 6 months, since these were modeled identically in GBD 2021.
ebf_df<-merge(sev_ebf[,-"age_group_id"], ebf_prev[age_group_id==3,], by=c("location_id", "year_id", "sex_id"))
ebf_df<- merge(ebf_df, sdi[,c("location_id", "year_id", "sdi", "lower_sdi", "upper_sdi")], by=c("location_id", "year_id"))
cor(ebf_df$prevalence, ebf_df$sev)
plot(ebf_df$sev, ebf_df$prevalence)
#change prevalence of EBF to prevalence of non-ebf by subtracting from 1
ebf_df[,prevalence_non:=1-prevalence]
ebf_df[,sev_ebf:= 1-sev ]
#change the REI name
ebf_df[,rei_name:="Exclusive breastfeeding"]
cor(ebf_df$prevalence, ebf_df$sev_ebf)

#plot for ebf
non_ebf_plot <- ggplot(ebf_df, aes(x=sev, y=prevalence_non, color=region_name)) +
  geom_point()+
  labs(title = "non-EBF proportion vs. non-EBF SEV (1990-2021)") +
  xlab("SEV") +
  ylab("Proportion non-EBF")+
  theme(legend.title = element_text(size = 4), legend.text = element_text(size = 4)) +
  geom_abline()+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  labs(subtitle=paste0("Correlation= ", formatC(cor(ebf_df$prevalence_non, ebf_df$sev),digits=3, format="f")) )

ebf_plot <- ggplot(ebf_df, aes(x=sev_ebf, y=prevalence, color=region_name)) +
  geom_point()+
  labs(title = "EBF proportion vs. (1 - non-EBF SEV) (1990-2021)") +
  xlab("1-SEV") +
  ylab("Proportion")+
  theme(legend.title = element_text(size = 4), legend.text = element_text(size = 4)) +
  geom_abline()+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  labs(subtitle=paste0("Correlation= ", formatC(cor(ebf_df$prevalence, ebf_df$sev_ebf),digits=3, format="f")) )



#wasting
#Prevalence estimates from GBD2021- include each of the neonatal agegroups
estimates1990_2021<-get_draws("modelable_entity_id", 10558, year_id=c(1990:2021),
                              source="epi", release_id = 9,
                              age_group_id = c(2,3,4,5), sex_id = c(1, 2, 3), location_id = locs)

#previously used the aggregated agegroup 1 (under 5s)
estimates1990_2021 <- melt(estimates1990_2021, id.vars = id.vars)
estimates1990_2021 <- prepare_df_format(df = estimates1990_2021, loc.met[level==3])
wast_prev <- merge(estimates1990_2021, hierarchy[, .(location_id, super_region_name, region_name)], all.x=T, by = "location_id")

wast_prev <- wast_prev[location_id%in%loc.met$location_id, c("location_id", "year_id", "age_group_id","sex_id","super_region_name","region_name","val", "lower","upper")]
setnames(wast_prev, old=c("val","lower", "upper"), new=c("prevalence","lower_prev", "upper_prev"))

populations <- get_population(age_group_id = c(388,389), location_set_id = 91, year_id = c(1990:2021), release_id = 9,location_id = locs, sex_id = c(1,2,3))

#population-weight SEVs from agegroup 388 and 389 to get a sev for agegroup 4 (post-neonatal)
sev_wasting_4<- sev_wasting_45[age_group_id%in%c(388,389),]
sev_wasting_4 <- merge(sev_wasting_4, populations[,-"run_id"], by= c("age_group_id", "location_id",  "year_id", "sex_id"))
sev_wasting_4[age_group_id>380, sev:= weighted.mean(val, population), by = c("location_id",  "measure_id",   "metric_id", 
                                                                              "rei_id", "sex_id", "year_id", "expected",
                                                                              "location_name", "location_type", "measure_name",  "metric_name",  
                                                                              "rei", "rei_name", "sex" )]
sev_wasting_4[age_group_id>380, lower_sev:= weighted.mean(lower, population), by = c("location_id",  "measure_id",   "metric_id", 
                                                                                      "rei_id", "sex_id", "year_id", "expected",
                                                                                      "location_name", "location_type", "measure_name",  "metric_name",  
                                                                                      "rei", "rei_name", "sex" )]
sev_wasting_4[age_group_id>380, upper_sev:= weighted.mean(upper, population), by = c("location_id",  "measure_id",   "metric_id", 
                                                                                      "rei_id", "sex_id", "year_id", "expected",
                                                                                      "location_name", "location_type", "measure_name",  "metric_name",  
                                                                                      "rei", "rei_name", "sex" )]
#drop val, lower, and upper
sev_wasting_4<- sev_wasting_4[age_group_id==388, -c("val","lower", "upper","population")]
sev_wasting_4[, age_group_id:=4]
sev_wasting_4[, age_group_name:= "Post Neonatal"]
#merge sev_wasting_4 with wast_prev and then population-weight the prevalence to get prevalence for agegroup 4
wast_df4<- merge(sev_wasting_4, wast_prev, by=c("location_id", "year_id", "sex_id", "age_group_id"))

#fill the SEVs from age_group_id==4 into SEV for age_group_id==2 and 3
wast_df23<-merge(sev_wasting_4[, -c("age_group_id","age_group_name")], wast_prev[age_group_id%in%c(2,3)], by=c("location_id", "year_id", "sex_id"))

#population-weight SEVs from agegroup 238 and 34 to get a sev for agegroup 5 (1 to 4 years)
sev_wasting_5<- sev_wasting_45[age_group_id%in%c(238,34),]

#setnames(sev_wasting_5, old=c("val","lower", "upper"), new=c("sev","lower_sev", "upper_sev"))
populations <- get_population(age_group_id = c(238,34), location_set_id = 91, year_id = c(1990:2021), release_id = 9, location_id = locs, sex_id = c(1,2,3))

#population-weight SEVs from agegroup 388 and 389 to get a sev for agegroup 4 (post-neonatal)
sev_wasting_5 <- merge(sev_wasting_5, populations[,-"run_id"], by= c("age_group_id", "location_id",  "year_id", "sex_id"))
sev_wasting_5[, sev:= weighted.mean(val, population), by = c("location_id",  "measure_id",   "metric_id", 
                                                              "rei_id", "sex_id", "year_id", "expected",
                                                              "location_name", "location_type", "measure_name",  "metric_name",  
                                                              "rei", "rei_name", "sex" )]
sev_wasting_5[, lower_sev:= weighted.mean(lower, population), by = c("location_id",  "measure_id",   "metric_id", 
                                                                      "rei_id", "sex_id", "year_id", "expected",
                                                                      "location_name", "location_type", "measure_name",  "metric_name",  
                                                                      "rei", "rei_name", "sex" )]
sev_wasting_5[, upper_sev:= weighted.mean(upper, population), by = c("location_id",  "measure_id",   "metric_id", 
                                                                      "rei_id", "sex_id", "year_id", "expected",
                                                                      "location_name", "location_type", "measure_name",  "metric_name",  
                                                                      "rei", "rei_name", "sex" )]
#select just one agegroup to sub in for age_group 5
sev_wasting_5<- sev_wasting_5[age_group_id==34, -c("val","lower", "upper", "population")]
sev_wasting_5[,age_group_id:=5]
sev_wasting_5[, age_group_name:= "1 to 4"]
#merge sev_wasting_4 with wast_prev and then population-weight the prevalence to get prevalence for agegroup 4
wast_df5<- merge(sev_wasting_5, wast_prev, by=c("location_id", "year_id", "sex_id", "age_group_id"))

#rbind together wast_df and wast_df_23
wast_df <- rbind(wast_df23, wast_df4, wast_df5, fill=TRUE)
wast_df[age_group_id == 2, age_group_name := "Early Neonatal"]
wast_df[age_group_id == 3, age_group_name := "Late Neonatal"]

wast_df<- merge(wast_df, sdi[,c("location_id", "year_id", "sdi", "lower_sdi", "upper_sdi")], by=c("location_id", "year_id"))
#correlation for agegroup2 prevalence with agegroup4 SEV
cor(wast_df[sex_id==2 & age_group_id==2]$prevalence, wast_df[sex_id==2 & age_group_id==2]$sev)
cor(wast_df[year_id==1990 & sex_id ==3]$prevalence, wast_df[year_id==1990 & sex_id==3]$sev)
cor(wast_df[year_id==2010 & sex_id ==3]$prevalence, wast_df[year_id==2010 & sex_id==3]$sev)
cor(wast_df[year_id==2019 & sex_id ==3]$prevalence, wast_df[year_id==2019 & sex_id==3]$sev)

plot(wast_df$sev, wast_df$prevalence)
#plot for wasting
wast_plot <- ggplot(wast_df[sex_id==3 & age_group_id==5], aes(x=sev, y=prevalence, color=region_name)) +
  geom_point()+
  labs(title = "wasting prevalence vs. SEV (1990- 2021), 1 to 4 years") +
  xlab("SEV") +
  ylab("Prevalence")+
  theme(legend.title = element_text(size = 4), legend.text = element_text(size = 4)) +
  geom_abline()+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  labs(subtitle=paste0("Correlation= ", formatC(cor(wast_df[sex_id==3 & age_group_id==5]$prevalence, wast_df[sex_id==3 & age_group_id==5]$sev),digits=3, format="f")) )


#MR-BRT models- prep the SEV dataframe to go into MR-BRT
#merge the different indicator DFs, but not anemia
sev_df<-rbind(stunt_df, lbw_df, ebf_df[,c(1:28)], wast_df)
# use offset before logit
sev_df[lower_prev<1e-7, lower_prev := 1e-7]
sev_df[upper_prev>(1-1e-7), upper_prev := (1-1e-7)]
#create logit_prevalence
sev_df[,logit_prev:= logit(prevalence)]
#this produces se (logit p) from lower and upper bounds 
sev_df[, se:= (logit(upper_prev) - logit(lower_prev))/3.92]


#prep the cgf dataframe to go into MR-BRT (until the other dfs are brought up to speed)
cgf_df<-rbind(stunt_df, wast_df)
# use offset before logit
cgf_df[lower_prev<1e-7, lower_prev := 1e-7]
cgf_df[upper_prev>(1-1e-7), upper_prev := (1-1e-7)]
#create logit_prevalence
cgf_df[,logit_prev:= logit(prevalence)]
#this produces se (logit p) from lower and upper bounds 
cgf_df[, se:= (logit(upper_prev) - logit(lower_prev))/3.92]


#prep the overweight dataframe to go into MR-BRT
# use offset before logit
overweight_df[lower_prev<1e-7, lower_prev := 1e-7]
overweight_df[upper_prev>(1-1e-7), upper_prev := (1-1e-7)]
#create logit_prevalence
overweight_df[,logit_prev:= logit(prevalence)]
#this produces se (logit p) from lower and upper bounds 
overweight_df[, se:= (logit(upper_prev) - logit(lower_prev))/3.92]
overweight_df[, age_group_id:=34]
overweight_df[age_group_id == 34, age_group_name := "2 to 4"]
overweight_df<- overweight_df[, c("location_id", "year_id",     "location_name","sex_id", "age_group_id", "measure_id",   "metric_id",   
                                "rei_id", "age_group_name", "expected", "location_type","measure_name", "metric_name",  "rei",     
                                "rei_name", "sex", "sev_kids", "upper_sev_kids", "lower_sev_kids", "super_region_name", "region_name", 
                                "prevalence",   "lower_prev",   "upper_prev",   "sev_adults", "lower_sev_adults", "upper_sev_adults", "sdi",     
                                "lower_sdi",    "upper_sdi",    "logit_prev",   "se" )]


#prep the anemia dataframe to go into MR-BRT
#create logit_prevalence
anem_df[,logit_prev:= logit(prevalence)]
#this produces se (logit p) from lower and upper bounds 
anem_df[, se:= (logit(upper_prev) - logit(lower_prev))/3.92]


saveRDS(sev_df, file.path(output_dir, "SEV_prev2021.RDS"))
sev_df <- readRDS(file.path(output_dir, "SEV_prev2021.RDS"))
saveRDS(overweight_df, file.path(output_dir, "overweight_prev2021.RDS"))
overweight_df <- readRDS(file.path(output_dir, "overweight_prev2021.RDS"))
saveRDS(anem_df, file.path(output_dir, "Hb_prev2021.RDS"))
anem_df <- readRDS(file.path(output_dir, "Hb_prev2021.RDS"))
saveRDS(cgf_df, file.path(output_dir, "cgf_prev2021.RDS"))
cgf_df <- readRDS(file.path(output_dir, "cgf_prev2021.RDS"))

#creating the agegroup-specific anemia dataframe:
estimates1990_2021<- get_draws("modelable_entity_id", 10507, year_id=c(1990:2021), location_id= locs, source="epi", release_id = 9, age_group_id = c(8:14), sex_id = 2)
estimates1990_2021 <- melt(estimates1990_2021, id.vars = id.vars)
estimates1990_2021 <- prepare_df_format(df = estimates1990_2021, loc.met[level==3])
anem_prev <- merge(estimates1990_2021, hierarchy[, .(location_id, super_region_name, region_name)], all.x=T, by = "location_id")

anem_prev <- anem_prev[location_id%in%loc.met$location_id & sex_id==2, c("location_id", "year_id", "age_group_id","super_region_name","region_name","val", "lower","upper")]

sev_irondef = get_outputs("rei", rei_id=95, measure_id=29, metric_id=3, location_id = locs, age_group_id = c(8:14, 24), sex_id = 2, release_id= 9, compare_version_id= 7999, year_id=c(1990:2021))

setnames(anem_prev, old=c("val","lower", "upper"), new=c("prevalence","lower_prev", "upper_prev"))
setnames(sev_irondef, old=c("val","lower", "upper"), new=c("sev","lower_sev", "upper_sev"), skip_absent=TRUE)
anem_df_all<-merge(sev_irondef, anem_prev, by=c("location_id", "year_id", "age_group_id"))
anem_df_all<- merge(anem_df_all, sdi[,c("location_id", "year_id", "sdi", "lower_sdi", "upper_sdi")], by=c("location_id", "year_id"))
anem_df_all[, logit_prev:= logit(prevalence)]
anem_df_all[, se:= (logit(upper_prev) - logit(lower_prev))/3.92]
saveRDS(anem_df_all, file.path(output_dir, "anem_prev2021.RDS"))
anem_df_all<- readRDS(file.path(output_dir, "anem_prev2021.RDS"))


#by indicator, run MR-BRT and cascading splines:
output_dir1<-"FILEPATH"
fit_global <- function(model_data, r) {
  #subset to only contain both sexes and split into training and testing dataframes
  #anemia will be run separately (sex_id==2)
  train_df<- model_data[year_id<2015 & rei_id==r,]
  test_df <- model_data[year_id>=2015 & rei_id==r,-"logit_prev"]
  
  #looking at sev distribution
  hist(model_data[rei_id==r,]$sev)
  #training data
  dat_loc <- MRData()
  dat_loc$load_df(
    data = train_df,
    col_obs = "logit_prev", col_obs_se = "se",
    col_covs = list("sdi", "sev"), col_study_id = "location_id"
  )
  #prediction data
  dat_pred <- MRData()
  dat_pred$load_df(
    data = test_df,
    col_covs=list("sdi", "sev"), col_study_id = "location_id"
  )

    #"standard MRBRT model
  stand_mod <- MRBRT(
    data = dat_loc,
    cov_models = list(
      LinearCovModel("intercept", use_re = TRUE),
      LinearCovModel("sdi", use_re = FALSE),
      LinearCovModel("sev", use_re = FALSE, 
                     use_spline = TRUE,
                     spline_knots = array(seq(0, 1, by = 0.25)), # this tells is to put a spline every 0.25
                     spline_degree = 3L,
                     spline_knots_type = 'frequency',# this specifies to put the knots literally along the domain. Other option is frequency for it to be data density driven
                     spline_r_linear = TRUE,
                     spline_l_linear = FALSE
      )
    ))
  stand_mod$fit_model(inner_print_level = 0L, inner_max_iter = 1000L,  outer_max_iter= 500L)
  estimated_beta <- stand_mod$summary()[[1]][1, "sdi"]
  stand_mod$cov_models[[which(stand_mod$cov_names == "sdi")]]$prior_beta_uniform <- matrix(rep(estimated_beta, 2), ncol = 1)
  #make predictions from global model
  pred_logit <- stand_mod$predict(data = dat_pred, predict_for_study= TRUE, sort_by_data_id= TRUE)
  test_df$pred_logit <- pred_logit
  test_df[,pred_prev:= inv.logit(pred_logit)]
  
  plot(test_df$pred_prev, test_df$prevalence, main= paste0(unique(test_df$rei_name), " predicted vs. observed prevalence in 2015-2021"), 
       sub= paste0("RMSE= ", rmse(test_df$prevalence, test_df$pred_prev)), xlab = "Predicted prevalence", ylab= "Observed prevalence")
  abline(0,1)
  
  print(paste0(unique(test_df$rei_name), " RMSE= ", rmse(test_df$prevalence, test_df$pred_prev)))
  
  return(stand_mod)
  
}

#need new global fit for stunting and wasting, to predict prevalence for all the under 5 age_groups
fit_global_cgf <- function(model_data, r, agid) {
  #subset to only contain both sexes and split into training and testing dataframes
  train_df<- model_data[year_id<2015 & rei_id==r & age_group_id==agid,]
  test_df <- model_data[year_id>=2015 & rei_id==r & age_group_id==agid,-"logit_prev"]
  
  #looking at sev distribution
  hist(model_data[rei_id==r & age_group_id==agid,]$sev)
  #training data
  dat_loc <- MRData()
  dat_loc$load_df(
    data = train_df,
    col_obs = "logit_prev", col_obs_se = "se",
    col_covs = list("sdi", "sev"), col_study_id = "location_id"
  )
  #prediction data
  dat_pred <- MRData()
  dat_pred$load_df(
    data = test_df,
    col_covs=list("sdi", "sev"), col_study_id = "location_id"
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
  
  stand_mod$fit_model(inner_print_level = 0L, inner_max_iter = 1000L,  outer_max_iter= 500L)
  estimated_beta <- stand_mod$summary()[[1]][1, "sdi"]
  stand_mod$cov_models[[which(stand_mod$cov_names == "sdi")]]$prior_beta_uniform <- matrix(rep(estimated_beta, 2), ncol = 1)
  #make predictions from global model
  pred_logit <- stand_mod$predict(data = dat_pred, predict_for_study= TRUE, sort_by_data_id= TRUE)
  test_df$pred_logit <- pred_logit
  test_df[,pred_prev:= inv.logit(pred_logit)]
  
  plot(test_df$pred_prev, test_df$prevalence, main= paste0(unique(test_df$rei_name), " in ", unique(test_df$age_group_name), " predicted vs. observed prevalence in 2015-2021"), 
       sub= paste0("RMSE= ", rmse(test_df$prevalence, test_df$pred_prev)), xlab = "Predicted prevalence", ylab= "Observed prevalence")
  abline(0,1)
  
  print(paste0(unique(test_df$rei_name), " in ", unique(test_df$age_group_name), " RMSE= ", rmse(test_df$prevalence, test_df$pred_prev)))
  
  return(stand_mod)
  
}


fit_global_overweight <- function(model_data) {
  #subset to only contain both sexes and split into training and testing dataframes
  train_df<- model_data[year_id<2015,]
  test_df <- model_data[year_id>=2015,-"logit_prev"]
  
  #looking at sev distribution
  hist(model_data$sev_adults)
  #training data
  dat_loc <- MRData()
  dat_loc$load_df(
    data = train_df,
    col_obs = "logit_prev", col_obs_se = "se",
    col_covs = list("sdi", "sev_kids", "year_id"), col_study_id = "location_id"
  )
  #prediction data
  dat_pred <- MRData()
  dat_pred$load_df(
    data = test_df,
    col_covs=list("sdi", "sev_kids", "year_id"), col_study_id = "location_id"
  )

  #"standard MRBRT model
  stand_mod <- MRBRT(
    data = dat_loc,
    cov_models = list(
      LinearCovModel("intercept", use_re = TRUE),
      LinearCovModel("sdi", use_re = FALSE),
      LinearCovModel("sev_kids", use_re = FALSE, 
                     use_spline = TRUE,
                     #fit spline to every quartile of the SEV distribution
                     spline_knots = array(seq(0, 1, by = 0.25)), # this tells is to put a spline every 0.25
                     spline_degree = 3L,
                     spline_knots_type = 'frequency',# this specifies to put the knots literally along the domain. Other option is frequency for it to be data density driven
                     spline_r_linear = TRUE,
                     spline_l_linear = FALSE
      )
    ))
  
  stand_mod$fit_model(inner_print_level = 5L, inner_max_iter = 500L,  outer_max_iter= 1000L)
  
  #make predictions from global model
  pred_logit <- stand_mod$predict(data = dat_pred, predict_for_study= TRUE, sort_by_data_id= TRUE)
  test_df$pred_logit <- pred_logit
  test_df[,pred_prev:= inv.logit(pred_logit)]
  
  plot(test_df$pred_prev, test_df$prevalence, main= paste0(unique(test_df$rei_name), " predicted vs. observed prevalence in 2015-2021"), 
       sub= paste0("RMSE= ", rmse(test_df$prevalence, test_df$pred_prev)), xlab = "Predicted prevalence", ylab= "Observed prevalence")
  abline(0,1)
  
  print(paste0(unique(test_df$rei_name), " RMSE= ", rmse(test_df$prevalence, test_df$pred_prev)))
  
  return(stand_mod)
  
}

fit_global_anemia <- function(model_data, agid) {
  #subset to only contain both sexes and split into training and testing dataframes
  train_df<- model_data[year_id<2015 & age_group_id== agid,]
  test_df <- model_data[year_id>=2015 & age_group_id== agid,-"logit_prev"]
  
  #training data
  dat_loc <- MRData()
  dat_loc$load_df(
    data = train_df,
    col_obs = "logit_prev", col_obs_se = "se",
    col_covs = list("sdi", "year_id", "sev"), col_study_id = "location_id"
  )
  #prediction data
  dat_pred <- MRData()
  dat_pred$load_df(
    data = test_df,
    col_covs=list("sdi", "year_id", "sev"), col_study_id = "location_id"
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
                     spline_l_linear = FALSE,
                     prior_spline_convexity = "concave"
       )
    ))
  
  stand_mod$fit_model(inner_print_level = 0L, inner_max_iter = 500L,  outer_max_iter= 1000L)
  #do not pass prior_beta_uniform for SDI for overweight because the relationship may be changing over time
   estimated_beta <- stand_mod$summary()[[1]][1, "sdi"]
   stand_mod$cov_models[[which(stand_mod$cov_names == "sdi")]]$prior_beta_uniform <- matrix(rep(estimated_beta, 2), ncol = 1)
   #make predictions from global model
  pred_logit <- stand_mod$predict(data = dat_pred, predict_for_study= TRUE, sort_by_data_id= TRUE)
  test_df$pred_logit <- pred_logit
  test_df[,pred_prev:= inv.logit(pred_logit)]
  
  plot(test_df$pred_prev, test_df$prevalence, main= paste0("Anemia in WRA: predicted vs. observed prevalence in 2015-2021"), 
       sub= paste0("RMSE= ", rmse(test_df$prevalence, test_df$pred_prev)), xlab = "Predicted prevalence", ylab= "Observed prevalence")
  abline(0,1)
  
  print(paste0(unique(test_df$rei_name), " in ", unique(test_df$age_group_name), " RMSE= ", rmse(test_df$prevalence, test_df$pred_prev)))
  
  return(stand_mod)
}

#cascading spline functions:

#cascading spline for EBF and LBW
fit_cascade <- function(mod_global, model_data, output_dir, r) {

    train_df<- model_data[year_id<2015 & rei_id==r,]

  model_label_tmp <- paste0("cascade_prevalence", unique(model_data[rei_id==r,]$rei_id))
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

#cascading spline for child growth failure risks
fit_cascade_cgf <- function(mod_global, model_data, output_dir, r, agid) {
  
  train_df<- model_data[year_id<2015 & rei_id==r & age_group_id== agid,]
  
  model_label_tmp <- paste0("cascade_prevalence", unique(train_df[rei_id==r,]$rei_id), "_", unique(train_df[rei_id==r,]$age_group_id))
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

#cascading spline for BMI and anemia
fit_cascade_2ind <- function(mod_global, model_data, output_dir, agid) {
 
    train_df<- model_data[year_id<2015 & age_group_id== agid,]

#BMI
   if(unique(train_df$rei_id)==371 ){
     name<-"high_bmi_SEV_kids"
    model_label_tmp <- paste0("cascade_prevalence", unique(model_data$rei_id), name)
    theta_list<- c(1,3)
   }
    
#anemia/iron deficiency    
   if(unique(train_df$rei_id)==95 ){
    name<-"iron_SEV"
    model_label_tmp <- paste0("cascade_prevalence", unique(train_df$rei_id), "_", unique(train_df$age_group_id))
    theta_list<- c(5,3)
  }

  thetas <- theta_list
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



#fit CGF age-specific fit_global models and cascading splines
for(r in c(241,240)){

for(agid in unique(cgf_df$age_group_id)){
  if(r==241)    {
    assign(paste0("stunting", agid), fit_global_cgf(cgf_df[sex_id==3], 241, agid))
    assign(paste0("stunt_cascade", agid), fit_cascade_cgf(get(paste0("stunting", agid)),cgf_df[sex_id==3],output_dir1, 241, agid ))
  }
  if(r==240)    {
    assign(paste0("wasting", agid), fit_global_cgf(cgf_df[sex_id==3], 240, agid))
    assign(paste0("waste_cascade", agid), fit_cascade_cgf(get(paste0("wasting", agid)),cgf_df[sex_id==3],output_dir1, 240, agid ))
  }
}
}

#fit global models for LBW, EBF, overweight, anemia
lbw <- fit_global(sev_df[sex_id==3], 335)
ebf <- fit_global(sev_df[sex_id==3],136)
overweight <- fit_global_overweight(overweight_df[sex_id==3 & !is.na(sev_kids)])
anemia <- fit_global_anemia(anem_df, 24)
for(agid in unique(anem_df_all$age_group_id)){
    assign(paste0("anemia", agid), fit_global_anemia(anem_df_all, agid))
    assign(paste0("anemia_cascade", agid), fit_cascade_2ind(get(paste0("anemia", agid)),anem_df_all,output_dir1, agid ))
}

lbw_cascade <- fit_cascade(lbw, sev_df[sex_id==3], output_dir1, 335)
ebf_cascade <- fit_cascade(ebf,sev_df[sex_id==3], output_dir1, 136)
overweight_cascade <- fit_cascade_2ind(overweight,overweight_df[sex_id==3 & !is.na(sev_kids)], output_dir1, 34)
anemia_cascade<- fit_cascade_2ind(anemia,anem_df,output_dir1)
