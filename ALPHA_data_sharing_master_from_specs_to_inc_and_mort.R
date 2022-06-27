##########################################################################
### DATA SHARING MASTER R SCRIPT
##########################################################################

# PREPARES DATA FROM ALPHA RESIDENCY AND HIV TEST SPECS TO PRODUCE
# DATASETS ON HIV MORTALITY AND HIV INCIDENCE FOR SHARING 
# VIA DATAFIRST
# OCTOBER 2021
# 
# Before running this R script, need to ensure all data and R scripts (supplied as a zip file)
# in the correct directories.
# Also, if not using post-neg times from metadata, need to ensure you have 
# post-neg times in a folder structure 
# {alphapath}/ALPHA/Estimates_Incidence/Post_negative_times/

if(!require(tidyverse)) install.packages("tidyverse")
if(!require(haven)) install.packages("haven")
if(!require(survival)) install.packages("survival")
if(!require(labelled)) install.packages("labelled")
if(!require(rlang)) install.packages("rlang")
library(tidyverse); library(haven); library(survival); library(labelled); library(rlang)

# +=+=+=  THINGS YOU NEED TO SET BEFORE RUNNING THE R SCRIPT +=+=+=+=
# set sitename - you must use the ALPHA sitename, as used in the filenames within the zip file
sitename <- "Karonga"
# set the path to the drive/directory where you have placed the ALPHA folder with data in the 
# appropriate the sub-folders
alphapath <- "K:/Katie/Do_file_translation_to_R"
# set the location of the zip file 

# prefix for create_hiv_status_detail script
prefix <- ""

##################################################
##################################################
source(paste0(alphapath,"/Prepare_data/prepare_residency_summary_dates.R"))
source(paste0(alphapath,"/Prepare_data/Get_HIV_data_from_hiv_tests_check_residency_make_wide_for_merge.R"))

source(paste0(alphapath,"/Analysis/Make_analysis_file_incidence_midpoint.R"))
source(paste0(alphapath,"/Analysis/Make_analysis_file_incidence_mi.R"))

source(paste0(alphapath,"/Analysis/Make_analysis_file_mortality_by_HIV_status.R"))
