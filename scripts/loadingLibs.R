# alorenzetti 202012

# description ####
# this script will load all
# the libs required to perform
# the differential expression analysis
# and load essential variables

# loading libs ####
# biocmanager is required for loading and installing packages
if(!requ ire("BiocManager")){install.packages("BiocManager"); library("BiocManager")}

# pacman is a nice package manager; make it easier to load and install packages
if(!require("pacman")){install.packages("pacman"); library("pacman")}

# list of required packages
requiredPacks = c("tidyverse",
                  "pheatmap",
                  "RColorBrewer",
                  "genefilter",
                  "DT",
                  "ggplot2",
                  "plotly",
                  "Rqc",
                  "rtracklayer",
                  "QuasR",
                  "GenomicAlignments",
                  "Rsamtools",
                  "GenomicFeatures",
                  "BiocParallel",
                  "DESeq2",
                  "gage",
                  "openxlsx",
                  "magrittr",
                  "UniProt.ws")

# loading required cran packages
p_load(char=requiredPacks)
