# alorenzetti 202012

# description ####
# this script will get and wrangle annotation
# and functional data for genes

# setting up annotation ####
# loading annotation file
gffFile = file.path("../misc/Ccrescentus.gff")
annot = rtracklayer::import(gffFile)

# filtering genes
genes = subset(annot, type == "gene")
genes$Dbxref = genes$Dbxref[lapply(genes$Dbxref, grepl, pattern = "GeneID:")] %>%
  sub(pattern = "GeneID:", replacement = "") %>%
  as.numeric()

# filtering CDS to obtain protein products
CDS = subset(annot, type == "CDS")
CDS$Dbxref = CDS$Dbxref[lapply(CDS$Dbxref, grepl, pattern = "GeneID:")] %>%
  sub(pattern = "GeneID:", replacement = "") %>%
  as.character()

CDS = CDS[CDS$Dbxref %>% duplicated == FALSE,]
CDS = as_tibble(CDS)
CDS = CDS %>%
  dplyr::select(Dbxref, product)

names(genes) = paste(genes$ID, genes$locus_tag, genes$Dbxref, genes$Name, sep = "|")

# getting kegg info ####
if(file.exists("data/keggSetCcsDf.RData")){
  load("data/keggSetCcsDf.RData")
  
} else {
  keggSetCcs = kegg.gsets(species = "ccs", id.type = "entrez")
  keggSetCcs = keggSetCcs$kg.sets
  
  keggSetCcsDf = keggSetCcs %>%
    unlist %>%
    tibble::enframe()
  
  keggSetCcsDf$name = sub(keggSetCcsDf$name, pattern = "[0-9]{1,}$", replacement = "")
  colnames(keggSetCcsDf) = c("KEGGpathway", "entrezid")
  
  keggSetCcsDf %<>%
    dplyr::group_by(entrezid) %>%
    dplyr::summarise(KEGGpathway = base::paste(KEGGpathway, collapse = "; "))
  
  save(keggSetCcsDf, file="data/keggSetCcsDf.RData")
}

# getting uniprot info ####
if(file.exists("data/uniprotCcs.RData")){
  load("data/uniprotCcs.RData")
  
} else {
  uniprotCcs = UniProt.ws(taxId=565050)
  columns = c("UNIPROTKB", "GO")
  uniprotCcs = UniProt.ws::select(uniprotCcs,
                                  keys=CDS$Dbxref,
                                  columns = columns,
                                  keytype = "ENTREZ_GENE")
  
  uniprotCcs = uniprotCcs[uniprotCcs$ENTREZ_GENE %>% duplicated == FALSE,]
  
  save(uniprotCcs, file="data/uniprotCcs.RData")
}

# getting cog info ####
# cog only exists for ccrescentus cb15
# we are working with na1000
# ive got corresponding cb15 gis to na1000 locus_tag
# from ortholuge db
# to infer cog of na1000 based on cb15
if(file.exists("data/cogCcs.RData")){
  load("data/cogCcs.RData")
} else {
  # reading ortholuge
  gi = read_delim("http://www.pathogenomics.sfu.ca/ortholugedb/paired/download/plain?strain1=Caulobacter%20crescentus%20CB15&strain2=Caulobacter%20crescentus%20NA1000", delim="\t")
  
  gi = gi %>% 
    dplyr::select(gicb15 = "GI (Strain 1)",
                  ltna1000 = "Locus Tag (Strain 2)") %>% 
    dplyr::mutate(gicb15 = as.character(gicb15),
                  ltna1000 = as.character(ltna1000)) %>% 
    dplyr::group_by(ltna1000) %>% 
    filter(row_number()==1)
  
  # reading cog
  cogNames = read_delim("ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab", delim = "\t")
  cogClasses = read_delim("ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/fun2003-2014.tab", delim = "\t")
  # it is going to throw a warning
  # since file has a comma at EOL
  cogData = read_csv("ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cog2003-2014.csv",
                     col_names = c("domain_id", "genome_name", "protein_id",
                                   "protein_length", "domain_start", "domain_end",
                                   "COG_id", "membership_class"))
  
  # cogData could have more than one COG for a GI
  # sometimes there is even a repetitive COG for the same GI
  # due to distinct domain regions
  cogData = cogData %>% 
    filter(genome_name == "Caulobacter_crescentus_CB15_uid57891") %>% 
    dplyr::select(domain_id, COG_id) %>% 
    dplyr::distinct()
  
  # parsing cogNames
  # be careful, it parses max four func codes per entry
  cogNames$func = cogNames$func %>%
    sub("(.)(.)", "\\1,\\2", .) %>%
    sub("(.),(.)(.)", "\\1,\\2,\\3", .) %>%
    sub("(.),(.),(.)(.)", "\\1,\\2,\\3,\\4", .)
  cogNames = cogNames %>%
    separate_rows(., func, sep = ",")
  
  cogFunct = left_join(cogNames, cogClasses, by = c("func" = "# Code")) %>% 
    dplyr::select("# COG", "func", "Name") %>% 
    dplyr::group_by(`# COG`) %>% 
    summarise(funcCode = paste0(func, collapse = "|"),
              COG = paste0(Name, collapse = "|"))
  
  cogFinal = left_join(cogData, cogFunct, by = c("COG_id" = "# COG")) %>% 
    mutate(COG = COG,
           COGcode = funcCode) %>% 
    dplyr::select(gi = domain_id,
                  COG_id, 
                  COGcode,
                  COG) %>% 
    group_by(gi, COGcode, COG) %>% 
    summarise(COG_id = paste0(COG_id, collapse = ";")) %>% 
    ungroup() %>% 
    group_by(gi) %>% 
    summarise(COGcode = paste0(COGcode, collapse = ";"),
              COG = paste0(COG, collapse = ";"),
              COG_id = paste0(COG_id, collapse = ";")) %>% 
    mutate(gi = as.character(gi))
  
  cogFinal = left_join(cogFinal,
                       cogNames %>% dplyr::select(-func),
                       by = c("COG_id" = "# COG")) %>% 
    rename(name = "COGproduct")
  
  cogCcs = left_join(cogFinal, gi, by = c("gi"  = "gicb15")) %>% 
    dplyr::mutate(locus_tag = ltna1000) %>% 
    dplyr::select(locus_tag, COGproduct, COG, COGcode, COG_id)
  
  save(cogCcs, file="data/cogCcs.RData")
}

# declaring functional categorization function ####
functCat = function(sigTable){
  funcat = sigTable[,c(1:4,6,10)] %>%
    dplyr::mutate(upORdown = case_when(log2FoldChange <= -log2fcthreshold ~ "DownRegulated",
                                       log2FoldChange >= log2fcthreshold ~ "Upregulated")) %>% 
    dplyr::left_join(CDS, by=c("entrezid" = "Dbxref")) %>% 
    dplyr::left_join(keggSetCcsDf, by="entrezid") %>% 
    dplyr::left_join(uniprotCcs, by=c("entrezid" = "ENTREZ_GENE")) %>% 
    dplyr::left_join(cogCcs, by="locus_tag") %>% 
    dplyr::arrange(log2FoldChange)
  
  funcat$KEGGpathway = funcat$KEGGpathway %>%
    replace_na("Undefined")
  funcat$UNIPROTKB = funcat$UNIPROTKB %>%
    replace_na("Undefined")
  funcat$GO = funcat$GO %>%
    replace_na("Undefined")
  funcat$COG = funcat$COG %>%
    replace_na("Undefined")
  funcat$COGcode = funcat$COGcode %>%
    replace_na("Undefined")
  funcat$COG_id = funcat$COG_id %>%
    replace_na("Undefined")
  funcat$COGproduct = funcat$COGproduct %>% 
    replace_na("Undefined")
  funcat$product = funcat$product %>% 
    replace_na("Undefined")
  
#  funcat$UNIPROTKB = paste0("<a href='https://www.uniprot.org/uniprot/",funcat$UNIPROTKB,"'>", funcat$UNIPROTKB,"</a>")
  
  return(funcat)
}

# creating a functional category table for all genes
funCat = genes %>%
  as_tibble() %>%
  dplyr::select(locus_tag, Name, Dbxref) %>%
  dplyr::mutate(Dbxref = as.character(Dbxref)) %>% 
  dplyr::left_join(CDS, by=c("Dbxref" = "Dbxref")) %>% 
  dplyr::left_join(keggSetCcsDf, by=c("Dbxref" = "entrezid")) %>% 
  dplyr::left_join(uniprotCcs, by=c("Dbxref" = "ENTREZ_GENE")) %>% 
  dplyr::left_join(cogCcs, by="locus_tag") %>% 
  dplyr::rename(entrezID = "Dbxref") %>% 
  dplyr::mutate(across(.cols = everything(),
                       .fns = ~ case_when(is.na(.x) ~ "Undefined",
                                          TRUE ~ as.character(.x)))) %>% 
  dplyr::distinct()

# saving tables to store funCat object
write.table(x = funCat,
            file = paste0("results/", "proteinFunctionalCategorization", ".tsv"),
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t",
            dec = ",")