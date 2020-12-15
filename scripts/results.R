# alorenzetti 202012

# description ####
# this script will declare
# functions to make plots
# and report the results

# exploratory analysis ####
# plotting heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# plotting principal component analysis
plotPCA(rld, intgroup = c("strain", "condition"))

# heatmap of gene clustering (based on rlog distance)
# the 100 genes with highest variance between samples
topVarGenes = head(order(rowVars(assay(rldNonBlind)),decreasing=TRUE),100)
mat = assay(rldNonBlind)[ topVarGenes, ]
mat = mat - rowMeans(mat)
df = as.data.frame(colData(rldNonBlind)[,c("condition","strain")])
pheatmap(mat, annotation_col=df, width=7, height = 7.5)

# defining interactive volcano plot function ####
volcaPlot = function(allTable, sigTable){
  # Volcano plot with plotly
  ax = list(
    title = "Log2(FC)",
    zeroline = TRUE,
    showline = TRUE,
    showticklabels = TRUE,
    showgrid = TRUE
  )
  
  ay = list(
    title = "-Log10(Adjusted p-value)",
    zeroline = TRUE,
    showline = TRUE,
    showticklabels = TRUE,
    showgrid = TRUE
  )
  
  # if there are no DE genes we don't have a way to plot them
  # so this chunk will handle it
  if(dim(sigTable)[1] == 0){
    plot_ly(allTable,
            x = allTable$log2FoldChange,
            y = -log10(allTable$padj),
            text=allTable$geneName,
            type="scatter",
            mode="markers",
            color="#4E79A7") %>%
      add_trace(x = log2fcthreshold, mode="lines", text=paste("Log_2(FC) =",log2fcthreshold),
                line = list(color = "black", width = 1)) %>%
      add_trace(x = -log2fcthreshold, mode="lines", text=paste("Log_2(FC) =",-log2fcthreshold),
                line = list(color = "black", width = 1)) %>%
      add_trace(y = -log10(padjthreshold), mode="lines", text="Adjusted p-value = 0.01",
                line = list(color = "black", width = 1)) %>%
      add_trace(y = -log10(padjthreshold), mode="lines", text=paste("Adjusted p-value =", padjthreshold),
                line = list(color = "black", width = 1)) %>%
      layout(showlegend=FALSE, xaxis = ax, yaxis = ay)
  }else{
    plot_ly(allTable,
            x = allTable$log2FoldChange,
            y = -log10(allTable$padj),
            text=allTable$geneName,
            type="scatter",
            mode="markers") %>%
      add_markers(x = sigTable$log2FoldChange,
                  y = -log10(sigTable$padj),
                  text=sigTable$geneName,
                  mode="markers", type="scatter",
                  marker=list(color="#E15759")) %>%
      add_trace(x = log2fcthreshold, mode="lines", text=paste("Log_2(FC) =",log2fcthreshold),
                line = list(color = "black", width = 1)) %>%
      add_trace(x = -log2fcthreshold, mode="lines", text=paste("Log_2(FC) =",-log2fcthreshold),
                line = list(color = "black", width = 1)) %>%
      add_trace(y = -log10(padjthreshold), mode="lines", text="Adjusted p-value = 0.01",
                line = list(color = "black", width = 1)) %>%
      add_trace(y = -log10(padjthreshold), mode="lines", text=paste("Adjusted p-value =", padjthreshold),
                line = list(color = "black", width = 1)) %>%
      layout(showlegend=FALSE, xaxis = ax, yaxis = ay)
  }
}

# reporting results ####
# defining function to report results
generateResults = function(resdf){
  
  # converting dds to tibble
  resdf = resdf %>% 
    as_tibble(rownames = "gffid|locus_tag|entrezid|geneName") %>% 
    separate(col = `gffid|locus_tag|entrezid|geneName`, sep = "\\|", into = c("gffid", "locus_tag", "entrezid", "geneName"))
  
  # getting significant genes
  resdfsig = resdf %>% 
    dplyr::filter(abs(log2FoldChange) >= log2fcthreshold & padj < padjthreshold)
  
  # getting upregulated ones
  resdfup = resdf %>% 
    dplyr::filter(log2FoldChange >= log2fcthreshold) %>% 
    dplyr::arrange(desc(log2FoldChange)) %>%
    mutate_if(is.numeric, signif, digits = 4) %>% 
    datatable(rownames = FALSE, options = list(pageLength = 5))
  
  # getting downregulated ones
  resdfdown = resdf %>%
    dplyr::filter(log2FoldChange <= -log2fcthreshold) %>% 
    dplyr::arrange(log2FoldChange) %>%
    mutate_if(is.numeric, signif, digits = 4) %>% 
    datatable(rownames = FALSE, options = list(pageLength = 5))
  
  # getting a tibble containing all genes
  resdfall = resdf %>% 
    mutate_if(is.numeric, signif, digits = 4) %>% 
    datatable(rownames = FALSE, options = list(pageLength = 5))
  
  # generating interactive volcano plots
  resdfvolcano = volcaPlot(resdf, resdfsig)
  
  # performing functional categorization for significant genes
  resdffuncat = functCat(resdfsig) %>% 
    mutate_if(is.numeric, signif, digits = 4) %>% 
    datatable(escape=F, rownames = FALSE, options = list(pageLength = 5))
  
  # creating and organizing list to store
  # previous objects
  reslist = list()
  reslist[["upregulated"]] = resdfup
  reslist[["downregulated"]] = resdfdown
  reslist[["all"]] = resdfall
  reslist[["funcat"]] = resdffuncat
  reslist[["volcano"]] = resdfvolcano
  
  return(reslist)
}

# generating results ####
contrasts = c(
  "RhlB30C_vs_NA100030C",
  "RhlB10C_vs_NA100010C",
  "NA100010C_vs_NA100030C"
)

# getting objects
finRes = list()
for(i in contrasts){
  finRes[[i]] = generateResults(results[[i]])
}

# writing tables ####
contrasts = c(
  "RhlB30C_vs_NA100030C",
  "RhlB10C_vs_NA100010C",
  "NA100010C_vs_NA100030C"
)

for(i in contrasts){
  allGenes = i
  sigGenes = paste0(i, "Sig")
  funcatGenes = paste0(i, "Sig_funcat")
  
  write.table(as_tibble(get0(allGenes)), paste0("data/", allGenes, ".tsv"), col.names = T, row.names = T, quote = F, sep = "\t", dec = ",")
  write.xlsx(as_tibble(get0(allGenes)), paste0("data/", allGenes, ".xlsx"))
  
  write.table(as_tibble(get0(sigGenes)), paste0("data/", sigGenes, ".tsv"), col.names = T, row.names = T, quote = F, sep = "\t", dec = ",")
  write.xlsx(as_tibble(get0(sigGenes)), paste0("data/", sigGenes, ".xlsx"))
  
  write.table(as_tibble(get0(funcatGenes) %>%
                          mutate(UNIPROTKB = sub(UNIPROTKB, pattern="^.*/(.*)'>.*$", replacement = "\\1"))), paste0("data/", funcatGenes, ".tsv"), col.names = T, row.names = T, quote = F, sep = "\t", dec = ",")
  write.xlsx(as_tibble(get0(funcatGenes) %>%
                         mutate(UNIPROTKB = sub(UNIPROTKB, pattern="^.*/(.*)'>.*$", replacement = "\\1"))), paste0("data/", funcatGenes, ".xlsx"))
}

# writing whole count matrix ####
write.table(as_tibble(assay(se), rownames = "gffid|locus_tag|entrezid|geneName"), paste0("data/", "countMatrix", ".tsv"), col.names = T, row.names = T, quote = F, sep = "\t", dec = ",")
write.xlsx(as_tibble(assay(se), rownames = "gffid|locus_tag|entrezid|geneName"), paste0("data/", "countMatrix", ".xlsx"))

# summary
sigtables = paste0(contrasts, "Sig")

summaryTable=NULL

for(i in sigtables){
  up = dim(subset(get0(i), log2FoldChange > 0))[1]
  down = dim(subset(get0(i), log2FoldChange < 0))[1]
  summaryTable=rbind(summaryTable, c(i,up,down))
}

summaryTable=as.data.frame(summaryTable)
colnames(summaryTable) = c("sigTable", "up", "down")