################################################################################
## DESCRIPTION: Launch script for array job to run models to predict prevalence of Global Nutrition Target (GNT) indicators from Summary Exposure Values and socio-demographic index
## Before running this script to launch models, follow instructions here: https://github.com/ihmeuw-msca/mrtoolR

##############################################################################################################################################
##############################################################################################################################################
########################################################### SET UP ###########################################################################
##############################################################################################################################################
##############################################################################################################################################
library(readr)
library(data.table)
reis <- c(95, 240, 241, 335, 371, 136)
ages <- c(2,3,4,5, 34, 164, 8:14)
sexes <- c(1, 2)
scenarios<- 0

# create every unique combination of location, age, sex, scenario
param_map <- expand.grid(rei_id = reis,
                         age_group_id = ages,
                         sex_id = sexes, 
                         scenario = scenarios)
param_map <- data.table(param_map)
#subset REIs to appropriate ages
param_map <- param_map[(rei_id ==136 & age_group_id==3) | (rei_id == 335 & age_group_id==164) | (rei_id == 95 & sex_id==2 & age_group_id%in%c(8:14) | (rei_id == 371 & age_group_id== 34) | ((rei_id == 241 | rei_id == 240) & age_group_id%in%c(2:5))) ,]
param_map<- param_map[order(rei_id, sex_id, age_group_id, scenario)]

param_map_filepath <- "FILEPATH"

write_csv(param_map, param_map_filepath)
##############################################################################################################################################
##############################################################################################################################################
############################################################ Q SUB ###########################################################################
##############################################################################################################################################
##############################################################################################################################################
## QSUB Command
job_name <- "GNT_draws"   # name of the job
thread_flag <- "-c 3" 
mem_flag <- "--mem=20G" 
runtime_flag <- "-t 1000" #minutes allowed to run
queue_flag <- "-p long.q" # long or all

throttle_flag <- as.character(24) # how many tasks are allowed to run at once. 800 is a good limit for smaller jobs. Bigger jobs will need smaller throttles
n_jobs <- paste0("1-", nrow(param_map), "%", throttle_flag) # this means you're running one task for every row of the param map you made. 
#n_jobs <- paste0("{1-", nrow(param_map), "%", throttle_flag, "}") # this means you're running one task for every row of the param map you made. 
prev_job <- "nojobholds" # never change this
next_script <- "FILEPATH" # filepath to the script you want to launch. parallel_gnt_draws.R
error_filepath <- paste0("-e FILEPATH", Sys.getenv("USER"), "/errors/%x.e%j") # where are errors going to be saved
output_filepath <- paste0("-o FILEPATH", Sys.getenv("USER"), "/output/%x.o%j") # where are outputs going to be saved


project_flag<- "-A proj_htwt"  # make sure this is one you have permissions for

qsub_command2 <- paste( "sbatch", "-J", job_name, project_flag, mem_flag, thread_flag, runtime_flag, queue_flag, "-a", n_jobs, error_filepath, output_filepath,  "FILEPATH/execRscript.sh -i FILEPATH/ihme_rstudio_3630.img -s", next_script, param_map_filepath )

system(qsub_command2) # this is the go button
