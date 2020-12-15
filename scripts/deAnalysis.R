# alorenzetti 202012

# description ####
# this script will prepare data
# generate the counts, and
# perform differential expression
# analysis

# preparing data ####
# providing info about data
samples = c("NA1_10C", "NA2_10C", "NA3_10C",
            "NA1_30C", "NA2_30C", "NA3_30C",
            "RhlB1_10C", "RhlB2_10C", "RhlB3_10C",
            "RhlB1_30C", "RhlB2_30C", "RhlB3_30C")

reps = rep(1:3, 4)
strains = c(rep("NA1000", 6), rep("RhlB", 6))
conditions = c(rep("10C", 3), rep("30C", 3),
               rep("10C", 3), rep("30C", 3))

# info datatable
info = data.frame(replicate = reps, strain = strains, condition = conditions)
rownames(info) = samples

# providing the path for bam files
filenamesPaired = paste0("../coverage_tmp/", rownames(info), "-paired", ".bam")
filenamesUnpaired = paste0("../coverage_tmp/", rownames(info), "-unpaired", ".bam")

# setting bam file list using Rsamtools
bamfilesPaired = BamFileList(filenamesPaired)
bamfilesUnpaired = BamFileList(filenamesUnpaired)


# creating count tables (SummarizedExperiment) ####
if(file.exists("data/sePaired.RData")){
  load("data/sePaired.RData")
} else{
  sePaired = summarizeOverlaps(features=genes, reads=bamfilesPaired,
                               mode="IntersectionNotEmpty", # this parameter deserves attention
                               singleEnd=FALSE,
                               ignore.strand=FALSE,
                               fragments=FALSE,
                               preprocess.reads=invertStrand)
  save(sePaired, file = "data/sePaired.RData")
}

if(file.exists("data/seUnpaired.RData")){
  load("data/seUnpaired.RData")
} else{
  seUnpaired = summarizeOverlaps(features=genes, reads=bamfilesUnpaired,
                                 mode="IntersectionNotEmpty", # this parameter deserves attention
                                 singleEnd=TRUE,
                                 ignore.strand=FALSE,
                                 preprocess.reads=invertStrand)
  save(seUnpaired, file = "data/seUnpaired.RData")
}

se = sePaired
assay(se) = assay(sePaired) + assay(seUnpaired)

# giving info to colData 
colData(se) = DataFrame(info)

# creating DESeq object for entire experiment ####
# to do exploratory analysis and differential expression analysis
dds = se
dds$group = factor(paste0(dds$strain, "_", dds$condition))
dds = DESeqDataSet(dds, design = ~ group)

# doing rlog transformation for distance and PCA
# before eliminating genes with zero counts
rld = rlog(dds)
rldNonBlind = rlog(dds, blind=F)

# remove genes with zero counts
dds = dds[rowSums(counts(dds)) > 1, ]
dds = DESeq(dds)
#resultsNames(dds)

# setting distance matrix for dds
sampleDists = dist(t(assay(rld)))
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste( rld$strain, rld$condition, sep="_" )
colnames(sampleDistMatrix) = NULL
colors = colorRampPalette(rev(brewer.pal(9, "Reds")) )(255)

# result tables for contrasts
results = list()
results[["RhlB30C_vs_NA100030C"]] = results(dds, contrast= c("group", "RhlB_30C", "NA1000_30C"), alpha = padjthreshold)
results[["RhlB10C_vs_NA100010C"]] = results(dds, contrast= c("group", "RhlB_10C", "NA1000_10C"), alpha = padjthreshold)
results[["NA100010C_vs_NA100030C"]] = results(dds, contrast= c("group", "NA1000_10C", "NA1000_30C"), alpha = padjthreshold)
