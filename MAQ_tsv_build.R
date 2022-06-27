#clean enviorment and load in libraries
rm(list = ls(all.names = TRUE))
library(tidyverse)
#set working directory
setwd("/Users/michael_pavia/BinaRena")
#dastool
a<-read.delim("dastool.map",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header = F,col.names=c("contig","dastool"))
#GC
b<-read.delim("GC_contigs.txt",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header = F,col.names=c("contig","GC"))
b<-separate(b,contig,c("contig","delete"),sep = " ", extra = "merge")
#maxbin
h<-read.delim("maxbin.map",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header = F,col.names=c("contig","maxbin"))
#metabat
i<-read.delim("metabat.map",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header = F,col.names=c("contig","metabat"))
#kraken taxa
j<-read.delim("rank_names.tsv",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
colnames(j)<-c("contig","domain","phylum","class","order","family","genus","species")
#depth
k<-read.delim("MAQ_depth.abund",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = " ",header =T)
k<-k[,-2]
colnames(k)<-c("contig","MAQ_000_10_MG","MAQ_000_20_MG","MAQ_050_10_MG","MAQ_050_20_MG","MAQ_100_10_MG","MAQ_100_20_MG")
k<-semi_join(k,b)#keep only observations in df1 that match in df2.
#kmer 4
o<-read.delim("kmer4.pca.tsv",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
colnames(o)<-c("contig","PC1_k4","PC2_k4","PC3_k4","PC4_k4","PC5_k4")
p<-read.delim("kmer4.tsne.tsv",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
colnames(p)<-c("contig","tSNE1_k4","tSNE2_k4")
q<-read.delim("kmer4.umap.tsv",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
colnames(q)<-c("contig","UMAP1_k4","UMAP2_k4")
o<-o%>%
  full_join(p)%>%
  full_join(q)
#kmer 5
r<-read.delim("kmer5.pca.tsv",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
colnames(r)<-c("contig","PC1_k5","PC2_k5","PC3_k5","PC4_k5","PC5_k5")
s<-read.delim("kmer5.tsne.tsv",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
colnames(s)<-c("contig","tSNE1_k5","tSNE2_k5")
t<-read.delim("kmer5.umap.tsv",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
colnames(t)<-c("contig","UMAP1_k5","UMAP2_k5")
r<-r %>%
  full_join(s) %>%
  full_join(t)
#kmer 6
u<-read.delim("kmer6.pca.tsv",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
colnames(u)<-c("contig","PC1_k6","PC2_k6","PC3_k6","PC4_k6","PC5_k6")
v<-read.delim("kmer6.tsne.tsv",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
colnames(v)<-c("contig","tSNE1_k6","tSNE2_k6")
w<-read.delim("kmer6.umap.tsv",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header =T)
colnames(w)<-c("contig","UMAP1_k6","UMAP2_k6")
u<-u %>%
  full_join(v) %>%
  full_join(w)
#size
x<-read.csv("MAQ_contig_length.csv")
x<-x %>%
  `colnames<-`(c("contig","length")) %>%
  filter(length>=2000)
#put it all together
l<-b%>%select(-delete)%>%full_join(a)%>%full_join(h)%>%full_join(i)%>%full_join(j)%>%full_join(k)%>%full_join(o)%>%full_join(r)%>%full_join(u)%>%full_join(x)
aa<-read.csv("/bin_Stats/MAQ_binning_metrics.csv")
aa<-aa[,c(1:4,49:53)]
aa[is.na(aa)]<-""
for (a in unique(aa$binning_plan)) {
  b<-aa[aa$binning_plan==a,]
  if (a=="dastool") {
    print("dastool")
    b<-b[,-2]
    colnames(b)<-c("dastool","completeness_das","redundancy_das","phylum_das","class_das","order_das","family_das","genus_das")
    c<-full_join(b,l)
  }
  if (a=="maxbin") {
    print("maxbin")
    b<-b[,-2]
    colnames(b)<-c("maxbin","completeness_max","redundancy_max","phylum_max","class_max","order_max","family_max","genus_max")
    d<-full_join(b,l)
  }
  if (a=="metabat") {
    print("metabat")
    b<-b[,-2]
    colnames(b)<-c("metabat","completeness_met","redundancy_met","phylum_met","class_met","order_met","family_met","genus_met")
    e<-full_join(b,l)
  }
}
z<-c%>%
  full_join(d)%>%
  full_join(e)
#add in nitrogen gene potential
a<-read.delim("ko00910.lst",quote = "",row.names = NULL,stringsAsFactors = FALSE,sep = "\t",header = T)
b<-read.csv("../MAQ_koscan_clean.csv")
c<-semi_join(b,a) %>% 
  separate(gene.name,c("k127","contig_r","drop"),sep = "_", extra = "drop") %>% 
  mutate(contig=paste0(k127,"_",contig_r)) %>%
  select(contig,KO)%>%
  unique()%>%
  group_by(contig)%>%
  summarize(N_KO = str_c(KO, collapse = ","))
d<-semi_join(b,a) %>% 
  separate(gene.name,c("k127","contig_r","drop"),sep = "_", extra = "drop") %>% 
  mutate(contig=paste0(k127,"_",contig_r)) %>%
  select(contig,KO)%>%
  group_by(contig)%>%
  summarise(N_count=n()) %>% 
  full_join(c) %>%
  full_join(z)
#add in sulfur gene potential
e<-read.csv("kosulfur.lst",header = F)
colnames(e)<-c("KO")
f<-semi_join(b,e)
f<-separate(f,gene.name,c("k127","contig","drop"),sep = "_", extra = "drop")
f$contig<-paste0(f$k127,"_",f$contig)
g<-f%>%select(contig,KO)%>%
  unique()%>%
  group_by(contig)%>%
  summarize(sulfur_KO = str_c(KO, collapse = ","))
h<-f %>%
  select(contig,KO) %>%
  group_by(contig) %>%
  summarise(count_sulfur=n()) %>% 
  full_join(g) %>%
  full_join(d)
#write binarena tsv
write.table(d,"MAQ_Binarena_Input.tsv",sep = "\t",quote = F,row.names = F)