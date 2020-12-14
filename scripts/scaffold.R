# alorenzetti 202012

# description ####
# this script is the scaffold
# running all the scripts
# required to perform
# the differential expression analysis

# loading required variables ####
# setting up number of threads
threads=8

# should we run quality control?
reportFlag="yes"

# sourcing ####
# loading libs
source("scripts/loadingLibs.R")

if(reportFlag == "yes"){
  source("scripts/qualityControl.R")
}
