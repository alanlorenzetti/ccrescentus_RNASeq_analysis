# alorenzetti 202012

# description ####
# this script is the scaffold
# running all the scripts
# required to perform
# the differential expression analysis

# loading required variables ####
# setting up number of threads
threads=8

# sourcing ####
# loading libs
source("loadingLibs.R")

# should we run quality control?
reportFlag="yes"
source("qualityControl.R")

