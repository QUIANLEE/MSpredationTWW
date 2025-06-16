library(MicrobiotaProcess)
library(phyloseq)
library(tidyverse)
library(RColorBrewer)
library(pacman)
library(magrittr)
library(reshape2)
library(ggsci)


setwd(dir = 'C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os')

rm(list = ls())

otu <- "featuretable-RelativeFrequency.qza"
rep <- "rep-seqs.qza"
tree <- "rooted-tree.qza"
tax <- "taxonomy.qza"
sample <- "sample-metadata-dnaos.tsv"
ps_dada2 <- import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,
                          mapfilename=sample,treeio=tree)

phytax <- get_taxadf(obj=ps_dada2, taxlevel=5)
p <- phyloseq::otu_table(phytax) %>% as.data.frame(header=NULL,check.names=T)
pt <- t(p) %>% as.data.frame()%>%mutate(sample=row.names(.))%>%select(sample,everything())%>% set_rownames(seq_len(nrow(.)))
metadata <- read.table(file="sample-metadata-dnaos.tsv",sep="\t",header=T,check.names=FALSE ,row.names=NULL) %>%select(sample,group)
mergedtable <- merge(pt,metadata,by="sample")%>%select(group,sample,everything())
mergedata1 <- mergedtable[order(mergedtable$group),]%>%t()%>% as.data.frame()
write.table (mergedata1,file ="level-5.txt", sep ="\t",quote = FALSE,col.names =FALSE)
mergedata2 <- t(mergedata1)%>% as.data.frame()%>%filter(group %in% c("OSc-28", "OSn-28"))%>%t()%>% as.data.frame()
mergedata2 <- rownames_to_column(mergedata2, var = "Row_Name")
mergedata2$Row_Name <- gsub("d5__", "", mergedata2$Row_Name)
mergedata2$Row_Name <- gsub("f__un_d2__", "Un_", mergedata2$Row_Name)
mergedata2$Row_Name <- gsub("f__un_d3__", "Un_", mergedata2$Row_Name)
mergedata2$Row_Name <- gsub("f__un_d4__", "Un_", mergedata2$Row_Name)
mergedata2<- column_to_rownames(mergedata2, var = "Row_Name")
mergedata2["group", ] <- gsub("OSc-28", "Day28(CHX)", mergedata2["group", ])
mergedata2["group", ] <- gsub("OSn-28", "Day28", mergedata2["group", ])
write.table (mergedata2,file ="Lefse-L5-cyc+-28.txt", sep ="\t",quote = FALSE,col.names =FALSE)

