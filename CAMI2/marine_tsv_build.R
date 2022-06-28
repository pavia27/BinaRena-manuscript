  rm(list = ls(all.names = TRUE))
  library(tidyverse)
  setwd("/Users/michaelpavia/Dropbox (ASU)/temp/")
  a<-read.delim("kmer4_marine.txt",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
  b<-read.delim("marine_depth.txt",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
  b<-b%>%filter(contigLen>=2000)
  c<-prcomp(a[,c(2:257)],scale. = T)
  d<-as.data.frame(c$x)
  e<-cbind(b,d)
  e<-e[,c(1:4,6,8,10,12,14,16,18,20,22,24:26)]
  f<-Sys.glob("*map")
  g<-read.delim("gsa_pooled_mapping_short.strains.binning.txt",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T) %>% filter(LENGTH>=2000)
  for (h in f) {
    i<-read.delim(h,quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =F,,col.names=c("contigName",h))
    e<-full_join(e,i)
  }
  e<-full_join(e,g)  
  write.table(e,"marine_Binarena_Input.tsv",sep = "\t",quote = F,row.names = F)
