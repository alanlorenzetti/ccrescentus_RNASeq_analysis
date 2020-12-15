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
# the 75 genes with highest variance between samples
topVarGenes = head(order(rowVars(assay(rldNonBlind)),decreasing=TRUE),75)
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
generateResults = function(resdf, name){
  
  # padj == 0 will be the mininum padj multiplied by 10E-2
  # otherwise it cannot be displayed on the volcano plot
  minpadj = resdf %>% as_tibble() %>% drop_na() %>% dplyr::select(padj) %>% filter(padj != 0) %>% min()
  fixpadj = minpadj * 10E-2
  
  # converting dds to tibble and fixing padj == 0
  resdf = resdf %>% 
    as_tibble(rownames = "gffid|locus_tag|entrezid|geneName") %>% 
    separate(col = `gffid|locus_tag|entrezid|geneName`, sep = "\\|", into = c("gffid", "locus_tag", "entrezid", "geneName")) %>% 
    mutate(padj = case_when(padj == 0 ~ fixpadj,
                            TRUE ~ as.double(padj)))
  
  # getting significant genes
  resdfsig = resdf %>% 
    dplyr::filter(abs(log2FoldChange) >= log2fcthreshold & padj < padjthreshold)
  
  # getting upregulated ones
  resdfup = resdf %>% 
    dplyr::filter(log2FoldChange >= log2fcthreshold & padj < padjthreshold) %>% 
    dplyr::arrange(desc(log2FoldChange)) %>%
    mutate_if(is.numeric, signif, digits = 4) #%>% 
#    datatable(rownames = FALSE, options = list(pageLength = 5))
  
  # getting downregulated ones
  resdfdown = resdf %>%
    dplyr::filter(log2FoldChange <= -log2fcthreshold & padj < padjthreshold) %>% 
    dplyr::arrange(log2FoldChange) %>%
    mutate_if(is.numeric, signif, digits = 4) #%>% 
#    datatable(rownames = FALSE, options = list(pageLength = 5))
  
  # getting a tibble containing all genes
  resdfall = resdf %>% 
    mutate_if(is.numeric, signif, digits = 4) #%>% 
#    datatable(rownames = FALSE, options = list(pageLength = 5))
  
  # generating interactive volcano plots
  resdfvolcano = volcaPlot(resdf, resdfsig)
  
  # performing functional categorization for significant genes
  resdffuncat = functCat(resdfsig) %>% 
    mutate_if(is.numeric, signif, digits = 4) %>% 
    dplyr::distinct()
#    datatable(escape=F, rownames = FALSE, options = list(pageLength = 5))
  
  # creating and organizing list to store
  # previous objects
  reslist = list()
  reslist[["sig"]] = resdfsig
  reslist[["upregulated"]] = resdfup
  reslist[["downregulated"]] = resdfdown
  reslist[["all"]] = resdfall
  reslist[["funcat"]] = resdffuncat
  reslist[["volcano"]] = resdfvolcano
  
  # writing a table containing all genes
  write.table(x = reslist[["all"]],
              file = paste0("results/", name, ".tsv"),
              col.names = T,
              row.names = F,
              quote = F,
              sep = "\t",
              dec = ",")
  write.xlsx(reslist[["all"]],
                       paste0("results/", name, ".xlsx"))
  
  # writing a table containing only significant genes
  write.table(x = reslist[["sig"]],
              file = paste0("results/", name, "_sig.tsv"),
              col.names = T,
              row.names = F,
              quote = F,
              sep = "\t",
              dec = ",")
  write.xlsx(reslist[["sig"]],
             paste0("results/", name, "_sig.xlsx"))
  
  # writing a table containing significant genes
  # with functional categorization
  write.table(x = reslist[["funcat"]],
              file = paste0("results/", name, "_funCat.tsv"),
              col.names = T,
              row.names = F,
              quote = F,
              sep = "\t",
              dec = ",")
  write.xlsx(reslist[["funcat"]],
             paste0("results/", name, "_funCat.xlsx"))
  
  # saving interactive volcano plot
  htmlwidgets::saveWidget(widget = reslist[["volcano"]],
                          file = paste0("results/", name, "_volcanoPlot.html"))
  
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
  finRes[[i]] = generateResults(results[[i]], i)
}

# writing whole count matrix ####
write.table(as_tibble(assay(se), rownames = "gffid|locus_tag|entrezid|geneName"),
            paste0("results/", "countMatrix", ".tsv"),
            col.names = T,
            row.names = T,
            quote = F,
            sep = "\t",
            dec = ",")
write.xlsx(as_tibble(assay(se), rownames = "gffid|locus_tag|entrezid|geneName"),
           paste0("results/", "countMatrix", ".xlsx"))

# summary ####
summaryTable=NULL
for(i in names(finRes)){
  up = dim(finRes[[i]]$upregulated)[1]
  down = dim(finRes[[i]]$downregulated)[1]
  summaryTable=rbind(summaryTable, c(i,up,down))
}

summaryTable=as.data.frame(summaryTable)
colnames(summaryTable) = c("sigTable", "up", "down")
