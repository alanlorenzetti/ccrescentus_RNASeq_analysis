# alorenzetti 20210304

# description ####
# this script will perform
# enrichment analysis inside
# upregulated and downregulated groups
# and will also plot figures

# starting processing ####
# creating and parsing objx ####
# we are going to search for enriched
# groups within NA1000 10C vs NA1000 30C 
# upregulated and downregulated set of genes
# I will have to arbitrarily keep only one
# category for each gene
# I am going to keep only the first one
enrichdf = list()
enrichdf[["NA100010C_vs_NA100030C"]][["downregulated"]] = left_join(x = finRes[["NA100010C_vs_NA100030C"]][["downregulated"]],
                                                                    y = funCat,
                                                                    by = "locus_tag") %>% 
  dplyr::select(locus_tag,
                geneName,
                log2FoldChange,
                padj,
                COG) %>% 
  dplyr::mutate(status = case_when(log2FoldChange >= 1 ~ "Upregulated",
                                   log2FoldChange <= -1 ~ "Downregulated"),
                COG = str_replace(COG, "\\|.*$", ""),
                COG = str_replace(COG, ";.*$", ""))

enrichdf[["NA100010C_vs_NA100030C"]][["upregulated"]] = left_join(x = finRes[["NA100010C_vs_NA100030C"]][["upregulated"]],
                                                                  y = funCat,
                                                                  by = "locus_tag") %>% 
  dplyr::select(locus_tag,
                geneName,
                log2FoldChange,
                padj,
                COG) %>% 
  dplyr::mutate(status = case_when(log2FoldChange >= 1 ~ "Upregulated",
                                   log2FoldChange <= -1 ~ "Downregulated"),
                COG = str_replace(COG, "\\|.*$", ""),
                COG = str_replace(COG, ";.*$", ""))

enrichdf[["NA100010C_vs_NA100030C"]][["all"]] = bind_rows(enrichdf[["NA100010C_vs_NA100030C"]][["downregulated"]],
                                                          enrichdf[["NA100010C_vs_NA100030C"]][["upregulated"]])

enrichdf[["NA100010C_vs_NA100030C"]][["allfig"]] = enrichdf[["NA100010C_vs_NA100030C"]][["all"]] %>% 
  group_by(COG, status) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  mutate(count = case_when(status == "Downregulated" ~ count * -1,
                           TRUE ~ as.numeric(count)))

enrichdf[["NA100010C_vs_NA100030C"]][["funcat"]] = funCat %>% 
  mutate(COG = str_replace(COG, "\\|.*$", ""))

# enrichment analysis ####
# regulation status as primary clusters
# hypergeometric enrichment test of
# COG var inside
# regulation status
vars = "COG"

# creating list to store results
enrich = list()
inputobj = enrichdf[["NA100010C_vs_NA100030C"]][["all"]]
funcatobj = enrichdf[["NA100010C_vs_NA100030C"]][["funcat"]]
for(j in inputobj$status %>% unique()){
  curRegGroup = inputobj %>%
    filter(status == j)
  for(k in vars){
    curRegGroupVec = curRegGroup %>% 
      dplyr::select(k) %>% 
      unlist(use.names = F)
    
    curRegGroupLvs = curRegGroupVec %>% 
      unique()
    
    for(l in curRegGroupLvs){
      wb = sum(curRegGroupVec == l)
      vecu = funcatobj %>% 
        dplyr::select(k) %>% 
        unlist(use.names = F)
      wu = sum(vecu == l)
      bu = sum(vecu != l)
      drawn = curRegGroupVec %>% length()
      
      pval = phyper(q= wb, m= wu, n= bu, k = drawn, lower.tail = F)
      
      tib = tibble(regRule = j,
                   criteria = k,
                   level = l,
                   pval = pval)
      enrich = bind_rows(enrich, tib)
    }
  }
}

# correcting pvalues using BH method
# filtering by pval
enrich$qval = p.adjust(enrich$pval, method = "BH")
enrich = enrich[enrich$qval < qthr,]

# adjusting dataset to include enrichment analysis
# on differential expression category plot 
enrichdf[["NA100010C_vs_NA100030C"]][["allfig"]] = enrichdf[["NA100010C_vs_NA100030C"]][["allfig"]] %>% 
  left_join(x = .,
            y = enrich,
            by = c("COG" = "level", "status" = "regRule")) %>%
  dplyr::mutate(enrichStatus = case_when(qval < 0.001 ~ "***",
                                  qval < 0.01 & qval >= 0.001 ~ "**",
                                  qval < 0.05 & qval >= 0.01 ~ "*",
                                  TRUE ~ NA_character_),
         COG = factor(COG, COG %>% unique() %>% rev()))

# plotting figure ####
enrichdf[["NA100010C_vs_NA100030C"]][["plot"]] = enrichdf[["NA100010C_vs_NA100030C"]][["allfig"]] %>% 
  ggplot(aes(x = count,
             fill = status,
             y = COG)) +
  geom_col() +
  geom_text(aes(label = enrichStatus,
                group = status),
            vjust = 0.75) +
  scale_fill_manual(name = "",
                    values = c("Upregulated" = "#E15759", 
                               "Downregulated" = "#4E79A7")) +
  scale_x_continuous(name = "Count",
                     limits = c(-100, 100),
                     labels = abs) +
  theme(text = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "bottom")

# saving
ggsave(plot = enrichdf[["NA100010C_vs_NA100030C"]][["plot"]],
       filename = "plots/enrichmentOfCategories.tiff",
       width = 7,
       height = 4,
       unit = "in",
       device = "tiff",
       dpi = 600)

ggsave(plot = enrichdf[["NA100010C_vs_NA100030C"]][["plot"]],
       filename = "plots/enrichmentOfCategories.pdf",
       width = 7,
       height = 4,
       unit = "in",
       device = "pdf",
       dpi = 600)

# volcano plot figure ####
theme_set(theme_bw())
volcanoplot10C_vs_30C = finRes$NA100010C_vs_NA100030C$all %>% 
  mutate(status = case_when(log2FoldChange <= -log2fcthreshold &
                              padj < padjthreshold ~ "Downregulated",
                            log2FoldChange >= log2fcthreshold &
                              padj < padjthreshold ~ "Upregulated",
                            TRUE ~ "nonDE")) %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(alpha = 0.75,
             show.legend = F) +
  geom_hline(yintercept = -log10(padjthreshold),
             size = 0.25,
             linetype = "dashed") +
  geom_vline(xintercept = c(-log2fcthreshold,log2fcthreshold),
             size = 0.25,
             linetype = "dashed") +
  scale_color_manual(values = c("Downregulated" = "#4E79A7",
                                "Upregulated" = "#E15759",
                                "nonDE" = "#79706E")) +
  scale_x_continuous(limits = c(-12,12),
                     breaks = seq(-12, 12, 2)) +
  xlab("Log<sub>2</sub>(FC)") +
  ylab("-Log<sub>10</sub>(Adjusted *p*-value)") +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        text = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

ggsave(filename = "plots/volcanoplot10C_vs_30C.tiff",
       dpi = 600,
       plot = volcanoplot10C_vs_30C,
       units = "in",
       device = "tiff",
       width = 8,
       height = 5)
