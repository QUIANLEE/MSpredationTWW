#################################用于分析lefse#######################################
library(MicrobiotaProcess)
library(phyloseq)
library(tidyverse)
library(RColorBrewer)
library(pacman)
library(magrittr)
library(reshape2)
library(ggsci)


setwd(dir = 'C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-cs')

rm(list = ls())

otu <- "featuretable-RelativeFrequency.qza"
rep <- "rep-seqs.qza"
tree <- "rooted-tree.qza"
tax <- "taxonomy.qza"
sample <- "sample-metadata-dnacs.tsv"
ps_dada2 <- import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,
                          mapfilename=sample,treeio=tree)



phytax <- get_taxadf(obj=ps_dada2, taxlevel=5)
p <- phyloseq::otu_table(phytax) %>% as.data.frame(header=NULL,check.names=T)
pt <- t(p) %>% as.data.frame()%>%mutate(sample=row.names(.))%>%select(sample,everything())%>% set_rownames(seq_len(nrow(.)))
metadata <- read.table(file="sample-metadata-dnacs.tsv",sep="\t",header=T,check.names=FALSE ,row.names=NULL) %>%select(sample,group)
mergedtable <- merge(pt,metadata,by="sample")%>%select(group,sample,everything())
mergedata1 <- mergedtable[order(mergedtable$group),]%>%t()%>% as.data.frame()
write.table (mergedata1,file ="level-5-ra.txt", sep ="\t",quote = FALSE,col.names =FALSE)
mergedata2 <- t(mergedata1)%>% as.data.frame()%>%filter(group %in% c("CSc-28", "CSn-28"))%>%t()%>% as.data.frame()
mergedata2 <- rownames_to_column(mergedata2, var = "Row_Name")
mergedata2$Row_Name <- gsub("d5__", "", mergedata2$Row_Name)
mergedata2$Row_Name <- gsub("f__un_d2__", "Un_", mergedata2$Row_Name)
mergedata2$Row_Name <- gsub("f__un_d3__", "Un_", mergedata2$Row_Name)
mergedata2$Row_Name <- gsub("f__un_d4__", "Un_", mergedata2$Row_Name)
mergedata2<- column_to_rownames(mergedata2, var = "Row_Name")
# Change "CSc-28" to "Day28(CHX)" and "CSn-28" to "Day28" in the 'group' row
mergedata2["group", ] <- gsub("CSc-28", "Day28(CHX)", mergedata2["group", ])
mergedata2["group", ] <- gsub("CSn-28", "Day28", mergedata2["group", ])
write.table (mergedata2,file ="lefse-L5-cyc+-28.txt", sep ="\t",quote = FALSE,col.names =FALSE)

