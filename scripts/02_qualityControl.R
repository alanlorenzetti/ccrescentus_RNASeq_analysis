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

# getting insights about raw data ####
# getting filenames for each lib 
files = list.files(path = "../raw", pattern = "*.fastq.gz", full.names = T)

# getting info about paired-end reads
pairs = rep(1:(length(files)/2),2)
pairs = sort(pairs)

# running QC
qual1 = rqcQA(files, pair = pairs, workers=threads, n = 1000000, sample=T)

# creating a directory to store reports
if(!dir.exists("reports")){system("mkdir reports")}

# reporting
rqcReport(qual1, outdir = "reports", file = "raw-quality-Rqc")

# getting insights about trimmed data ####
# getting filenames for paired libs 
files = list.files(path = "../trimmed", pattern = ".*-paired.*.fastq.gz", full.names = T)

# getting info about paired-end reads
pairs = rep(1:(length(files)/2),2)
pairs = sort(pairs)

# running QC
qual2 = rqcQA(files, pair = pairs, workers=threads, n = 1000000, sample=T)

# reporting
rqcReport(qual2, outdir = "reports", file = "trimmed-quality-Rqc")

# getting filenames for unpaired libs
files = list.files(path = "../trimmed", pattern = "*-unpaired.*.fastq.gz", full.names = T)

# running QC
qual3 = rqcQA(files, workers=threads, n = 1000000, sample=T)

# reporting
rqcReport(qual3, outdir = "reports", file = "trimmed-unpaired-quality-Rqc")
