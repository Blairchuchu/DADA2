#@Use normalized table to do Deseq2
#@Draw volcano plot and vanns diagram and compare differential ASVs between metnod of VST and RLE
#@
#@author:tzu-yu Chu
#@data:07.11.2021

#@-----------------------------------------------

# Deseq2 can't be zero,negtive and decimal

# To do differential seperately (3 groups)
  names1 = c("CO-2087","CO-2091","CO-2116",
             "CO-2348","CO-2355","CO-2370",        
             "CO-2371","CO-2388","CO-2395",        
             "CO-2398","CO-2399","CO-2414",         
             "CO-2423","CO-2426","CO-2429",         
             "CO-2438","CO-2495","CO-2519",         
             "CO-2543","CO-2560","CO-2602",         
             "CO-2653","CO-2677","CO-2684",         
             "CO-2699","CO-2828","CO-2916",         
             "CO-3178","CO-3195","CO-3215",        
             "CO-3231","CO-3242","CO-3255",         
             "CO-3271","CO-3275","CO-3289",         
             "CO-3457","CO-3470","CO-3473",         
             "CO-3479","CO-3493","CO-3506",
             "102-DAA0070-BL","102-DAA0075-BL","102-DAA0077-BL",
             "102-DAA0078-BL","102-DAA0080-BL","102-DAA0081-BL",    
             "102-DAA0084-BL","102-DAA0086-BL","102-DAA0087-BL",   
             "102-DAA0088-BL","102-DAA0089-BL","102-DAA0091-BL",   
             "102-DAA0095-BL","102-DAA0096-BL","102-DAA0097-BL",   
             "102-DAA0098-BL","102-DAA0099-BL","103-DAA-E045-BL", 
             "103-DAA-E047-BL","103-DAA-E048-BL","103-DAA-E050-BL",  
             "103-DAA-E051-BL","103-DAA-E053-BL","103-DAA-E057-BL",
             "103-DAA-E059-BL","103-DAA-E061-BL","103-DAA-E062-BL",
             "103-DAA-E063-BL","103-DAA-E066-BL","DAA-E010-BL",
             "DAA-E012-BL","DAA-E014-BL","DAA-E015-BL",
             "DAA-E016-BL","DAA-E020-BL","DAA-E021-BL",
             "DAA-E023-BL","DAA-E038-BL","DAA-E041-BL",
             "DAA-E043-BL","DAA0046-BL","DAA0060-BL")
  
  names2 = c("CO-2087","CO-2091","CO-2116",
             "CO-2348","CO-2355","CO-2370",        
             "CO-2371","CO-2388","CO-2395",        
             "CO-2398","CO-2399","CO-2414",         
             "CO-2423","CO-2426","CO-2429",         
             "CO-2438","CO-2495","CO-2519",         
             "CO-2543","CO-2560","CO-2602",         
             "CO-2653","CO-2677","CO-2684",         
             "CO-2699","CO-2828","CO-2916",         
             "CO-3178","CO-3195","CO-3215",        
             "CO-3231","CO-3242","CO-3255",         
             "CO-3271","CO-3275","CO-3289",         
             "CO-3457","CO-3470","CO-3473",         
             "CO-3479","CO-3493","CO-3506",
             "102-DAA-0075-SVR","102-DAA-0077-SVR","102-DAA-0078-SVR",
             "102-DAA0070-SVR","102-DAA0080-SVR","102-DAA0081-SVR", 
             "102-DAA0084-SVR","102-DAA0086-SVR","102-DAA0087-SVR", 
             "102-DAA0088-SVR","102-DAA0089-SVR","102-DAA0091-SVR", 
             "102-DAA0095-SVR","102-DAA0096-SVR","102-DAA0097-SVR", 
             "102-DAA0098-SVR","102-DAA0099-SVR","103-DAA-E045-SVR",
             "103-DAA-E047-SVR","103-DAA-E048-SVR","103-DAA-E050-SVR",
             "103-DAA-E051-SVR","103-DAA-E053-SVR","103-DAA-E057-SVR",
             "103-DAA-E059-SVR","103-DAA-E061SVR","103-DAA-E062-SVR",
             "103-DAA-E063-SVR","103-DAA-E066-SVR","DAA-E010-SVR",    
             "DAA-E012-SVR","DAA-E014-SVR","DAA-E015-SVR",    
             "DAA-E016-SVR","DAA-E020-SVR","DAA-E021-SVR",    
             "DAA-E023-SVR","DAA-E038-SVR","DAA-E041-SVR",  
             "DAA-E043-SVR","DAA0046-SVR","DAA0060-SVR")
  
  names3 = c("102-DAA0070-BL","102-DAA0075-BL","102-DAA0077-BL",
             "102-DAA0078-BL","102-DAA0080-BL","102-DAA0081-BL",    
             "102-DAA0084-BL","102-DAA0086-BL","102-DAA0087-BL",   
             "102-DAA0088-BL","102-DAA0089-BL","102-DAA0091-BL",   
             "102-DAA0095-BL","102-DAA0096-BL","102-DAA0097-BL",   
             "102-DAA0098-BL","102-DAA0099-BL","103-DAA-E045-BL", 
             "103-DAA-E047-BL","103-DAA-E048-BL","103-DAA-E050-BL",  
             "103-DAA-E051-BL","103-DAA-E053-BL","103-DAA-E057-BL",
             "103-DAA-E059-BL","103-DAA-E061-BL","103-DAA-E062-BL",
             "103-DAA-E063-BL","103-DAA-E066-BL","DAA-E010-BL",
             "DAA-E012-BL","DAA-E014-BL","DAA-E015-BL",
             "DAA-E016-BL","DAA-E020-BL","DAA-E021-BL",
             "DAA-E023-BL","DAA-E038-BL","DAA-E041-BL",
             "DAA-E043-BL","DAA0046-BL","DAA0060-BL",
             "102-DAA-0075-SVR","102-DAA-0077-SVR","102-DAA-0078-SVR",
             "102-DAA0070-SVR","102-DAA0080-SVR","102-DAA0081-SVR", 
             "102-DAA0084-SVR","102-DAA0086-SVR","102-DAA0087-SVR", 
             "102-DAA0088-SVR","102-DAA0089-SVR","102-DAA0091-SVR", 
             "102-DAA0095-SVR","102-DAA0096-SVR","102-DAA0097-SVR", 
             "102-DAA0098-SVR","102-DAA0099-SVR","103-DAA-E045-SVR",
             "103-DAA-E047-SVR","103-DAA-E048-SVR","103-DAA-E050-SVR",
             "103-DAA-E051-SVR","103-DAA-E053-SVR","103-DAA-E057-SVR",
             "103-DAA-E059-SVR","103-DAA-E061SVR","103-DAA-E062-SVR",
             "103-DAA-E063-SVR","103-DAA-E066-SVR","DAA-E010-SVR",    
             "DAA-E012-SVR","DAA-E014-SVR","DAA-E015-SVR",    
             "DAA-E016-SVR","DAA-E020-SVR","DAA-E021-SVR",    
             "DAA-E023-SVR","DAA-E038-SVR","DAA-E041-SVR",  
             "DAA-E043-SVR","DAA0046-SVR","DAA0060-SVR")
  

# vst

vst_tab_n2tn = vst_trans_count_tab[,names1]
vst_tab_n2svr = vst_trans_count_tab[,names2]
vst_tab_tnsvr = vst_trans_count_tab[,names3]

sample_info_tab_n2tn = sample_info_tab[names1,]
sample_info_tab_n2svr = sample_info_tab[names2,]
sample_info_tab_tn_svr = sample_info_tab[names3,]

# vst n2tn
vst_data = round(vst_tab_n2tn) #Rounding
vst_data = vst_data+ 1   # plus1 because table have to without zero to use Deseq2 
vst_data = data.frame(vst_data)
sample_info_tab_n2tn$Cohort = as.factor(sample_info_tab_n2tn$Cohort)
dds = DESeqDataSetFromMatrix(countData=vst_data, 
                             colData=sample_info_tab_n2tn, 
                             design=~Cohort)


dep = DESeq(dds)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
res_n2tn <- results(dep,contrast=c("Cohort", "non-HCV", "Treatment-naive"))
head(results(dep, tidy = TRUE)) # let's look at the results table
summary(res_n2tn) # summary of results
DEG_n2tn=as.data.frame(res_n2tn)
write.csv(x=DEG_n2tn,file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/res_vst_n2tn.csv") 

# Draw the volcano plot to view result of differential ASVs 
ggplot(DEG_n2tn) + 
geom_point(aes(x = log2FoldChange, y = -log10(pvalue))) + 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.5) + 
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.5) +
  theme(legend.position = "none", 
          plot.title = element_text(size = rel(1.5), hjust = 0.5), 
          axis.title = element_text(size = rel(1.25)))+
  ggtitle("vst_n2tn")

# Remove NA and filter de GENE
resOrdered=res_n2tn[order(res_n2tn$pvalue),] # Order
sum(res_n2tn$padj<0.1, na.rm=TRUE) # na.rm remove na
sum(res_n2tn$pvalue<0.05, na.rm=TRUE) #1970  
sum(res_n2tn$pvalue<0.05&(res_n2tn$log2FoldChange > 1 | res_n2tn$log2FoldChange < -1)) #880
sum(res_n2tn$pvalue<0.05&res_n2tn$log2FoldChange > 1) #265
sum(res_n2tn$pvalue<0.05&res_n2tn$log2FoldChange < -1) #615

# 1
# Crateria gets 1970 ASVs
resSig=subset(resOrdered, pvalue<0.05)  
DEG_n2tn=as.data.frame(resSig)
write.csv(x=resSig, file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/filter_p_vst_results.csv")
# 2
# After this crateria get 880 ASVs,( pvalue < 0.05, |log2FoldChange| > 1 )
resSig=subset(resOrdered, pvalue<0.05 & (log2FoldChange > 1 | log2FoldChange < -1))  
DEG_n2tn=as.data.frame(resSig)
diff_gene_DESeq2 <- row.names(DEG_n2tn) 
resdata <- merge(as.data.frame(res), as.data.frame(counts(dep, normalized=TRUE)),by = "row.names",sort=FALSE)
write.csv(x = resSig, file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/filter_P_fc_vst_results.csv")


# vst n2svr
vst_data = round(vst_tab_n2svr) #Rounding
vst_data = vst_data+ 1   #plus1 because table have to without zero to use Deseq2 
vst_data = data.frame(vst_data)
sample_info_tab_n2svr$Cohort = as.factor(sample_info_tab_n2svr$Cohort)
dds = DESeqDataSetFromMatrix(countData=vst_data, 
                             colData=sample_info_tab_n2svr, 
                             design=~Cohort)


dep = DESeq(dds)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
res_n2svr <- results(dep,contrast=c("Cohort", "non-HCV", "SVR"))
head(results(dep, tidy = TRUE)) # let's look at the results table
summary(res_n2svr) # summary of results
DEG_n2svr=as.data.frame(res_n2svr)
write.csv(x=DEG1,file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/res_vst_n2svr.csv")

# Draw the volcano plot to view result of differential ASVs 
ggplot(DEG_n2svr) + 
geom_point(aes(x = log2FoldChange, y = -log10(pvalue))) + 
xlab("log2 fold change") + 
ylab("-log10 adjusted p-value") + 
geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.5) + 
geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.5) +
theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25)))+
ggtitle("vst_n2svr") 

# Remove NA and filter de GENE
resOrdered=res_n2svr[order(res_n2svr$pvalue),] # Order
sum(res_n2svr$padj<0.1, na.rm=TRUE) # na.rm remove na
sum(res_n2svr$pvalue<0.05, na.rm=TRUE)  #2059
sum(res_n2svr$pvalue<0.05&(res_n2svr$log2FoldChange > 1 | res_n2svr$log2FoldChange < -1)) #920
sum(res_n2svr$pvalue<0.05&res_n2svr$log2FoldChange > 1 ) #265

# 1
# Crateria gets 2059 ASVs
resSig=subset(resOrdered, padj<0.1)  
DEG_n2svr=as.data.frame(resSig)
write.csv(x=resSig, file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/filter_p_vst_results.csv")

# 2
# After this crateria get 920 ASVs,( pvalue < 0.05, |log2FoldChange| > 1 )
resSig=subset(resOrdered, pvalue<0.05 & (log2FoldChange > 1 | log2FoldChange < -1))  
DEG_n2svr=as.data.frame(resSig)
diff_gene_DESeq2 <- row.names(DEG_n2svr) 
resdata <- merge(as.data.frame(res), as.data.frame(counts(dep, normalized=TRUE)),by = "row.names",sort=FALSE)
write.csv(x = resSig, file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/filter_P_fc_vst_results.csv")


# vst tn2svr
vst_data = round(vst_tab_tnsvr) #Rounding
vst_data = vst_data+ 1   #plus1 because table have to without zero to use Deseq2 
vst_data = data.frame(vst_data)
sample_info_tab_tn_svr$Cohort = as.factor(sample_info_tab_tn_svr$Cohort)
dds = DESeqDataSetFromMatrix(countData=vst_data, 
                             colData=sample_info_tab_tn_svr, 
                             design=~Cohort)


dep = DESeq(dds)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
res_tn2svr <- results(dep,contrast=c("Cohort", "Treatment-naive", "SVR"))
head(results(dep, tidy = TRUE)) # let's look at the results table
summary(res_tn2svr) # summary of results
DEG_tn2svr=as.data.frame(res_tn2svr)
write.csv(x=DEG_tn2svr,file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/res_vst_tn2svr.csv") 

# Draw the volcano plot to view result of differential ASVs 
ggplot(DEG_tn2svr) + 
  geom_point(aes(x = log2FoldChange, y = -log10(pvalue))) + 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.5) + 
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.5) +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25)))+
  ggtitle("vst_tn2svr") 

# Remove NA and filter de GENE
resOrdered=res_tn2svr[order(res_tn2svr$pvalue),] # Order
sum(res_tn2svr$padj<0.1, na.rm=TRUE) # na.rm remove na
sum(res_tn2svr$pvalue<0.05, na.rm=TRUE)   #146
sum(res_tn2svr$pvalue<0.05&(res_tn2svr$log2FoldChange > 1 | res_tn2svr$log2FoldChange < -1)) #21
sum(res_tn2svr$pvalue<0.05&res_tn2svr$log2FoldChange > 1) #6
# 1
# Crateria gets 146 ASVs
resSig=subset(resOrdered, padj<0.1)  
DEG2=as.data.frame(resSig)
write.csv(x=resSig, file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/filter_p_vst_results.csv")

# 2
# After this crateria get 21 ASVs,( pvalue < 0.05, |log2FoldChange| > 1 )
resSig=subset(resOrdered, pvalue<0.05 & (log2FoldChange > 1 | log2FoldChange < -1))  
DEG_tn2svr=as.data.frame(resSig)
diff_gene_DESeq2 <- row.names(DEG_tn2svr) 
resdata <- merge(as.data.frame(res_tn2svr), as.data.frame(counts(dep, normalized=TRUE)),by = "row.names",sort=FALSE)
write.csv(x = resSig, file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/filter_P_fc_vst_results.csv")


# Draw vanns digram: here try 3 methods

library(VennDiagram)

#1*

#UP&DOWN
venn.diagram(
  x = list(n2tn=rownames(res_n2tn[which(res_n2tn$pvalue<=0.05&(res_n2tn$log2FoldChange > 1 | res_n2tn$log2FoldChange < -1)),]),
           n2svr=rownames(res_n2svr[which(res_n2svr$pvalue<=0.05&(res_n2svr$log2FoldChange > 1 | res_n2svr$log2FoldChange < -1)),]),
           tn2svr=rownames(res_tn2svr[which(res_tn2svr$pvalue<=0.05&(res_tn2svr$log2FoldChange > 1 | res_tn2svr$log2FoldChange < -1)),])
           
  ),
  category.names = c("n2tn" , "n2svr" , "tn2svr"),
  filename = 'C:/Users/user/Desktop/HCVandNOR/R_dada2/IMG_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#134262", '#00A4DE', '#98A7B5'),
  fill = c(alpha("#134262",0.6), alpha('#00A4DE',0.6), alpha('#98A7B5',0.6)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.5,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#134262", '#00A4DE', '#98A7B5'),
  rotation = 1)
dev.off()
#2
library("gplots")
pdf("C:/Users/user/Desktop/HCVandNOR/R_dada2/venn3.pdf")
venn(list(n2tn=rownames(res_n2tn[which(res_n2tn$pvalue<=0.05&(res_n2tn$log2FoldChange > 1 | res_n2tn$log2FoldChange < -1) ),]),
          n2svr=rownames(res_n2svr[which(res_n2svr$pvalue<=0.05&(res_n2svr$log2FoldChange > 1 | res_n2svr$log2FoldChange < -1)),]),
          tn2svr=rownames(res_tn2svr[which(res_tn2svr$pvalue<=0.05&(res_tn2svr$log2FoldChange > 1 | res_tn2svr$log2FoldChange < -1) ),])))
dev.off()
#3
library("ggVennDiagram")
a = 
  ggVennDiagram(list(n2tn=rownames(res_n2tn[which(res_n2tn$pvalue<=0.05&(res_n2tn$log2FoldChange > 1 | res_n2tn$log2FoldChange < -1) <= 0.05),]),
                     n2svr=rownames(res_n2svr[which(res_n2svr$pvalue<=0.05&(res_n2svr$log2FoldChange > 1 | res_n2svr$log2FoldChange < -1) <= 0.05),]),
                     tn2svr=rownames(res_tn2svr[which(res_tn2svr$pvalue<=0.05&(res_tn2svr$log2FoldChange > 1 | res_tn2svr$log2FoldChange < -1) <= 0.05),])), label_alpha = 0)
dev.off()

#UP
venn.diagram(
  x = list(n2tn=rownames(res_n2tn[which(res_n2tn$pvalue<=0.05&res_n2tn$log2FoldChange > 1 ),]),
           n2svr=rownames(res_n2svr[which(res_n2svr$pvalue<=0.05&res_n2svr$log2FoldChange > 1),]),
           tn2svr=rownames(res_tn2svr[which(res_tn2svr$pvalue<=0.05&res_tn2svr$log2FoldChange > 1),])
           
  ),
  category.names = c("n2tn" , "n2svr" , "tn2svr"),
  filename = 'C:/Users/user/Desktop/HCVandNOR/R_dada2/down_vst_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#134262", '#00A4DE', '#98A7B5'),
  fill = c(alpha("#134262",0.6), alpha('#00A4DE',0.6), alpha('#98A7B5',0.6)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.5,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#134262", '#00A4DE', '#98A7B5'),
  rotation = 1)
dev.off()

#DOWN
venn.diagram(
  x = list(n2tn=rownames(res_n2tn[which(res_n2tn$pvalue<=0.05&res_n2tn$log2FoldChange < -1 ),]),
           n2svr=rownames(res_n2svr[which(res_n2svr$pvalue<=0.05&res_n2svr$log2FoldChange < -1),]),
           tn2svr=rownames(res_tn2svr[which(res_tn2svr$pvalue<=0.05&res_tn2svr$log2FoldChange < -1),])
           
  ),
  category.names = c("n2tn" , "n2svr" , "tn2svr"),
  filename = 'C:/Users/user/Desktop/HCVandNOR/R_dada2/down_vst_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#134262", '#00A4DE', '#98A7B5'),
  fill = c(alpha("#134262",0.6), alpha('#00A4DE',0.6), alpha('#98A7B5',0.6)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.5,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#134262", '#00A4DE', '#98A7B5'),
  rotation = 1)
dev.off()


# RLE
rle_tab_n2tn = rle_trans_counts_tab[,names1]
rle_tab_n2svr = rle_trans_counts_tab[,names2]
rle_tab_tnsvr = rle_trans_counts_tab[,names3]

sample_info_tab_n2tn = sample_info_tab[names1,]
sample_info_tab_n2svr = sample_info_tab[names2,]
sample_info_tab_tn_svr = sample_info_tab[names3,]


# rle n2tn
rle_data = round(rle_tab_n2tn) #Rounding
rle_data = rle_data+ 1   #plus1 because table have to without zero to use Deseq2 
rle_data = data.frame(rle_data)
sample_info_tab_n2tn$Cohort = as.factor(sample_info_tab_n2tn$Cohort)
dds = DESeqDataSetFromMatrix(countData=rle_data, 
                             colData=sample_info_tab_n2tn, 
                             design=~Cohort)


dep = DESeq(dds)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
rle_res_n2tn <- results(dep,contrast=c("Cohort", "non-HCV", "Treatment-naive"))
head(results(dep, tidy = TRUE)) # let's look at the results table
summary(rle_res_n2tn) # summary of results
head(rle_res_n2tn)
DEG_n2tn=as.data.frame(rle_res_n2tn)
write.csv(x=DEG_n2tn,file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/res_vst_n2tn.csv") 

# Draw the volcano plot to view result of differential ASVs
ggplot(DEG_n2tn) + 
  geom_point(aes(x = log2FoldChange, y = -log10(pvalue))) + 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.5) + 
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.5) +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25)))+
  ggtitle("rle_n2tn")

# Remove NA and filter de GENE
resOrdered=rle_res_n2tn[order(rle_res_n2tn$pvalue),] # Order
sum(rle_res_n2tn$padj<0.1, na.rm=TRUE) # na.rm remove na
sum(rle_res_n2tn$pvalue<0.001, na.rm=TRUE) #4309  
sum(rle_res_n2tn$pvalue<0.001&(rle_res_n2tn$log2FoldChange > 3 | rle_res_n2tn$log2FoldChange < -3)) #2836
# 1
# Crateria gets 4309 ASVs
resSig=subset(resOrdered, pvalue<0.05)  
DEG_n2tn=as.data.frame(resSig)
write.csv(x=resSig, file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/filter_p_vst_results.csv")

# 2
# After this crateria get 2836 ASVs,( pvalue < 0.05, |log2FoldChange| > 1 )
resSig=subset(resOrdered, pvalue<0.05 & (log2FoldChange > 1 | log2FoldChange < -1))  
DEG_n2tn=as.data.frame(resSig)
diff_gene_DESeq2 <- row.names(DEG_n2tn) 
resdata <- merge(as.data.frame(res), as.data.frame(counts(dep, normalized=TRUE)),by = "row.names",sort=FALSE)
write.csv(x = resSig, file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/filter_P_fc_vst_results.csv")


# rle n2svr
rle_data = round(rle_tab_n2svr) #Rounding
rle_data = rle_data+ 1   #plus1 because table have to without zero to use Deseq2 
rle_data = data.frame(rle_data)
sample_info_tab_n2svr$Cohort = as.factor(sample_info_tab_n2svr$Cohort)
dds = DESeqDataSetFromMatrix(countData=rle_data, 
                             colData=sample_info_tab_n2svr, 
                             design=~Cohort)


dep = DESeq(dds)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
rle_res_n2svr <- results(dep,contrast=c("Cohort", "non-HCV", "SVR"))
head(results(dep, tidy = TRUE)) # let's look at the results table
summary(rle_res_n2svr) # summary of results
DEG_n2svr=as.data.frame(rle_res_n2svr)
write.csv(x=DEG1,file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/res_vst_n2svr.csv")

# Draw the volcano plot to view result of differential ASVs
ggplot(DEG_n2svr) + 
  geom_point(aes(x = log2FoldChange, y = -log10(pvalue))) + 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  geom_vline(xintercept = c(-3, 3), color = 'gray', size = 0.5) + 
  geom_hline(yintercept = -log(0.001, 10), color = 'gray', size = 0.5) +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25)))+
  ggtitle("rle_n2svr") 

# Remove NA and filter de GENE
resOrdered=rle_res_n2svr[order(rle_res_n2svr$pvalue),] # Order
sum(rle_res_n2svr$padj<0.1, na.rm=TRUE) # na.rm remove na
sum(rle_res_n2svr$pvalue<0.05, na.rm=TRUE)  #4486
sum(rle_res_n2svr$pvalue<0.05&(rle_res_n2svr$log2FoldChange > 1 | rle_res_n2svr$log2FoldChange < -1)) #2939

# 1
# Crateria gets 4486 ASVs
resSig=subset(resOrdered, padj<0.1)  
DEG_n2svr=as.data.frame(resSig)
write.csv(x=resSig, file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/filter_p_vst_results.csv")

# 2
# After this crateria get 2939 ASVs,( pvalue < 0.05, |log2FoldChange| > 1 )
resSig=subset(resOrdered, pvalue<0.05 & (log2FoldChange > 1 | log2FoldChange < -1))  
DEG_n2svr=as.data.frame(resSig)
diff_gene_DESeq2 <- row.names(DEG_n2svr) 
resdata <- merge(as.data.frame(rle_res_n2svr), as.data.frame(counts(dep, normalized=TRUE)),by = "row.names",sort=FALSE)
write.csv(x = resSig, file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/filter_P_fc_vst_results.csv")


# rle tn2svr
rle_data = round(rle_tab_tnsvr) #Rounding
rle_data = rle_data+ 1   #plus1 because table have to without zero to use Deseq2 
rle_data = data.frame(rle_data)
sample_info_tab_tn_svr$Cohort = as.factor(sample_info_tab_tn_svr$Cohort)
dds = DESeqDataSetFromMatrix(countData=rle_data, 
                             colData=sample_info_tab_tn_svr, 
                             design=~Cohort)


dep = DESeq(dds)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
rle_res_tn2svr <- results(dep,contrast=c("Cohort", "Treatment-naive", "SVR"))
head(results(dep, tidy = TRUE)) # let's look at the results table
summary(rle_res_tn2svr) # summary of results
DEG_tn2svr=as.data.frame(rle_res_tn2svr)
write.csv(x=DEG_tn2svr,file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/res_vst_tn2svr.csv") 


# Draw the volcano plot to view result of differential ASVs
ggplot(DEG_tn2svr) + 
  geom_point(aes(x = log2FoldChange, y = -log10(pvalue))) + 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.5) + 
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.5) +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25)))+
  ggtitle("rle_tn2svr") 

# Remove NA and filter de GENE
resOrdered=rle_res_tn2svr[order(rle_res_tn2svr$pvalue),] # Order
sum(rle_res_tn2svr$padj<0.1, na.rm=TRUE) # na.rm remove na
sum(rle_res_tn2svr$pvalue<0.05, na.rm=TRUE)   #1798
sum(rle_res_tn2svr$pvalue<0.05&(rle_res_tn2svr$log2FoldChange > 1 | rle_res_tn2svr$log2FoldChange < -1)) #939
# 1
# Crateria gets 1798 ASVs
resSig=subset(resOrdered, padj<0.1)  
DEG2=as.data.frame(resSig)
write.csv(x=resSig, file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/filter_p_vst_results.csv")

# 2
# After this crateria get 939 ASVs,( pvalue < 0.05, |log2FoldChange| > 1 )
resSig=subset(resOrdered, pvalue<0.05 & (log2FoldChange > 1 | log2FoldChange < -1))  
DEG_tn2svr=as.data.frame(resSig)
diff_gene_DESeq2 <- row.names(DEG_tn2svr) 
resdata <- merge(as.data.frame(rle_res_tn2svr), as.data.frame(counts(dep, normalized=TRUE)),by = "row.names",sort=FALSE)
write.csv(x = resSig, file = "C:/Users/user/Desktop/HCVandNOR/R_dada2/filter_P_fc_vst_results.csv")


# Draw vanns digram
library(VennDiagram)

# UP&DOWN
venn.diagram(
  x = list(n2tn=rownames(rle_res_n2tn[which(rle_res_n2tn$pvalue<=0.05&(rle_res_n2tn$log2FoldChange > 1 | rle_res_n2tn$log2FoldChange < -1)),]),
           n2svr=rownames(rle_res_n2svr[which(rle_res_n2svr$pvalue<=0.05&(rle_res_n2svr$log2FoldChange > 1 | rle_res_n2svr$log2FoldChange < -1)),]),
           tn2svr=rownames(rle_res_tn2svr[which(rle_res_tn2svr$pvalue<=0.05&(rle_res_tn2svr$log2FoldChange > 1 | rle_res_tn2svr$log2FoldChange < -1)),])
           
  ),
  category.names = c("n2tn" , "n2svr" , "tn2svr"),
  filename = 'C:/Users/user/Desktop/HCVandNOR/R_dada2/rle_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#134262", '#00A4DE', '#98A7B5'),
  fill = c(alpha("#134262",0.6), alpha('#00A4DE',0.6), alpha('#98A7B5',0.6)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.5,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#134262", '#00A4DE', '#98A7B5'),
  rotation = 1)
dev.off()

#UP
venn.diagram(
  x = list(n2tn=rownames(rle_res_n2tn[which(rle_res_n2tn$pvalue<=0.05&rle_res_n2tn$log2FoldChange > 1 ),]),
           n2svr=rownames(rle_res_n2svr[which(rle_res_n2svr$pvalue<=0.05&rle_res_n2svr$log2FoldChange > 1 ),]),
           tn2svr=rownames(rle_res_tn2svr[which(rle_res_tn2svr$pvalue<=0.05&rle_res_tn2svr$log2FoldChange > 1 ),])
           
  ),
  category.names = c("n2tn" , "n2svr" , "tn2svr"),
  filename = 'C:/Users/user/Desktop/HCVandNOR/R_dada2/up_rle_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#134262", '#00A4DE', '#98A7B5'),
  fill = c(alpha("#134262",0.6), alpha('#00A4DE',0.6), alpha('#98A7B5',0.6)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.5,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#134262", '#00A4DE', '#98A7B5'),
  rotation = 1)
dev.off()

# DOWN
venn.diagram(
  x = list(n2tn=rownames(rle_res_n2tn[which(rle_res_n2tn$pvalue<=0.05&rle_res_n2tn$log2FoldChange < -1 ),]),
           n2svr=rownames(rle_res_n2svr[which(rle_res_n2svr$pvalue<=0.05&rle_res_n2svr$log2FoldChange < -1 ),]),
           tn2svr=rownames(rle_res_tn2svr[which(rle_res_tn2svr$pvalue<=0.05&rle_res_tn2svr$log2FoldChange < -1 ),])
           
  ),
  category.names = c("n2tn" , "n2svr" , "tn2svr"),
  filename = 'C:/Users/user/Desktop/HCVandNOR/R_dada2/down_rle_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#134262", '#00A4DE', '#98A7B5'),
  fill = c(alpha("#134262",0.6), alpha('#00A4DE',0.6), alpha('#98A7B5',0.6)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.5,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#134262", '#00A4DE', '#98A7B5'),
  rotation = 1)
dev.off()


# compare VST& RLE
# n2tn
venn.diagram(
  x = list(vst = rownames(res_n2tn[which(res_n2tn$pvalue<=0.05&(res_n2tn$log2FoldChange > 1 | res_n2tn$log2FoldChange < -1)),]),
           rle = rownames(rle_res_n2tn[which(rle_res_n2tn$pvalue<=0.001&(rle_res_n2tn$log2FoldChange > 3 | rle_res_n2tn$log2FoldChange < -3)),]))
  ,
  category.names = c("VST" , "RLE"),
  filename = 'C:/Users/user/Desktop/HCVandNOR/R_dada2/comparison_venn_n2tn_motify.png',
  output = TRUE ,
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  col=c("#006699", "#FF6666"),
  fill = c(alpha("#006699",0.6), alpha("#FF6666",0.6)),
  fontfamily = "sans"
)

# n2svr
venn.diagram(
  x = list(vst = rownames(res_n2svr[which(res_n2svr$pvalue<=0.05&(res_n2svr$log2FoldChange > 1 | res_n2svr$log2FoldChange < -1)),]),
           rle = rownames(rle_res_n2svr[which(rle_res_n2svr$pvalue<=0.05&(rle_res_n2svr$log2FoldChange > 1 | rle_res_n2svr$log2FoldChange < -1)),]))
  ,
  category.names = c("VST" , "RLE"),
  filename = 'C:/Users/user/Desktop/HCVandNOR/R_dada2/comparison_venn_n2svr.png',
  output = TRUE ,
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  col=c("#006699", "#FF6666"),
  fill = c(alpha("#006699",0.6), alpha("#FF6666",0.6)),
  fontfamily = "sans"
)

# tn2svr
venn.diagram(
  x = list(vst = rownames(res_tn2svr[which(res_tn2svr$pvalue<=0.05&(res_tn2svr$log2FoldChange > 1 | res_tn2svr$log2FoldChange < -1)),]),
           rle = rownames(rle_res_tn2svr[which(rle_res_tn2svr$pvalue<=0.05&(rle_res_tn2svr$log2FoldChange > 1 | rle_res_tn2svr$log2FoldChange < -1)),]))
  ,
  category.names = c("VST" , "RLE"),
  filename = 'C:/Users/user/Desktop/HCVandNOR/R_dada2/comparison_venn_tn2svr.png',
  output = TRUE ,
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  col=c("#006699", "#FF6666"),
  fill = c(alpha("#006699",0.6), alpha("#FF6666",0.6)),
  fontfamily = "sans"
)