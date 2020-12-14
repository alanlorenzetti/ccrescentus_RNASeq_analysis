# alorenzetti 202012

# description ####
# this script will perform
# quality control of
# RNA-Seq reads

# setting up a workaround for Rqc ####
# Rqc needs a workround for the moment https://support.bioconductor.org/p/91401/
# replace internal function code to avoid the bug
.readFrequency <- function (chunk)
{
  tbl <- table(as.character(sread(chunk)))
  count <- as.integer(tbl)
  hash <- names(tbl)
  data.frame(hash, count, stringsAsFactors = FALSE)
}
assignInNamespace(".readFrequency", .readFrequency, "Rqc")