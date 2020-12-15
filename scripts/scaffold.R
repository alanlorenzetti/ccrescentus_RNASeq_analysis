# alorenzetti 202012

# description ####
# this script is the scaffold
# running all the scripts
# required to perform
# the differential expression analysis

# loading required variables ####
# setting up number of threads
threads=10

# registering multiple processors
register(MulticoreParam(workers = threads))

# should we run quality control?
reportFlag="no"

# DE analysis thresholds
# DESeq2 adjusted pval
padjthreshold = 0.01

# DeSeq2 log2FoldChange
log2fcthreshold = 2

# creating data directory
if(!dir.exists("data")){dir.create("data")}

# creating results directory
if(!dir.exists("results")){dir.create("results")}

# sourcing ####
# loading libs
source("scripts/loadingLibs.R")
if(reportFlag == "yes"){source("scripts/qualityControl.R")}
source("scripts/functionalCategorization.R")
source("scripts/deAnalysis.R")
source("scripts/results.R")

