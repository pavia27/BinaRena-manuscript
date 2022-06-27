#clean enviorment and load in libraries
rm(list = ls(all.names = TRUE))
library(tidyverse)
#set working directory
setwd("/Users/michael_pavia/BinaRena")
o<-list("00003","00010","00045","00076","00078","00101","00125","00147","00156","00160","00510","00521","00523","00528","00538","00540","00560","00678","00715","06128","06163","06165","06168","50070","50076","50395","80129","80142","80152")
  for (p in o) {
    #dastool
    a<-read.delim(paste0(p,".dastool.map"),quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header = F,col.names=c("contig","dastool"))
    #GC
    b<-read.delim(paste0(p,"_GC"),quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header = F,col.names=c("contig","GC"))
    b<-separate(b,contig,c("contig","delete"),sep = " ", extra = "merge")
    #essential hmm
    g<-read.csv(paste0(p,"_clean.csv"),header = F)
    g<-separate(g,V1,c("cell","merge"),sep = "_", extra = "drop")%>%mutate(contig=paste0(cell,"_",merge))%>%select(contig,V3)%>%unique()%>%group_by(contig)%>%summarize(essential_gene=str_c(V3, collapse = ","))
    #maxbin
    h<-read.delim(paste0(p,".maxbin.map"),quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header = F,col.names=c("contig","maxbin"))
    #metabat
    i<-read.delim(paste0(p,".metabat.map"),quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header = F,col.names=c("contig","metabat"))
    #taxa
    j<-read.delim(paste0(p,"_rank_names.tsv"),quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
    colnames(j)<-c("contig","domain","phylum","class","order","family","genus","species")
    #depth
    k<-read.delim(paste0(p,".metabat.abund"),quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
    k<-k[,-c(4,5)] %>% filter(contigLen>=2000)
    colnames(k)<-c("contig","contigLen","totalAvgDepth")
    #kmer 4
    x<-read.delim(paste0(p,".kmer4.pca.tsv"),quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
    colnames(x)<-c("contig","PC1_k4","PC2_k4","PC3_k4","PC4_k4","PC5_k4")
    y<-read.delim(paste0(p,".kmer4.tsne.tsv"),quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
    colnames(y)<-c("contig","tSNE1_k4","tSNE2_k4")
    q<-read.delim(paste0(p,".kmer4.umap.tsv"),quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
    colnames(q)<-c("contig","UMAP1_k4","UMAP2_k4")
    x<-x%>%full_join(y)%>%full_join(q)
    #kmer 5
    r<-read.delim(paste0(p,".kmer5.pca.tsv"),quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
    colnames(r)<-c("contig","PC1_k5","PC2_k5","PC3_k5","PC4_k5","PC5_k5")
    s<-read.delim(paste0(p,".kmer5.tsne.tsv"),quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
    colnames(s)<-c("contig","tSNE1_k5","tSNE2_k5")
    t<-read.delim(paste0(p,".kmer5.umap.tsv"),quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
    colnames(t)<-c("contig","UMAP1_k5","UMAP2_k5")
    r<-r%>%full_join(s)%>%full_join(t)
    #kmer 6
    u<-read.delim(paste0(p,".kmer6.pca.tsv"),quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
    colnames(u)<-c("contig","PC1_k6","PC2_k6","PC3_k6","PC4_k6","PC5_k6")
    v<-read.delim(paste0(p,".kmer6.tsne.tsv"),quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
    colnames(v)<-c("contig","tSNE1_k6","tSNE2_k6")
    w<-read.delim(paste0(p,".kmer6.umap.tsv"),quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
    colnames(w)<-c("contig","UMAP1_k6","UMAP2_k6")
    u<-u%>%full_join(v)%>%full_join(w)
    #final
    l<-b%>%select(-delete)%>%full_join(a)%>%full_join(g)%>%full_join(h)%>%full_join(i)%>%full_join(j)%>%full_join(k)%>%full_join(x)%>%full_join(r)%>%full_join(u)
    m<-l[,c(1,3,8,9)]
    m$remove<-rowSums(is.na(m))
    m<-m[!m$remove==3,]
    n<-semi_join(l,m) #keep only observations in df1 that match in df2.
    n[is.na(n)]<-""
    write.table(n, paste0(p,"_Binarena_Input.tsv"),sep = "\t",quote = F,row.names = F)
  }
