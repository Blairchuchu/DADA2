#@16s amplicon
#@ Summarize the result of ASV table
#@ Normalize:within ,between sample
#@ Phyloseq:alpha-,beta- diversity, taxonomy composition
#@ Visualize prevalence of ASV
#@author:tzu-yu Chu
#@data:07.03.2021

#@============================================================================
# Draw ASV result
library("magrittr") #for %>%
library("janitor")
library("tidyverse")
library("dplyr")

# Input
track = read.csv ("C:/Users/user/Desktop/HCVandNOR/R_dada2/merged.out_track.csv")[-c(127,128,129),]
rownames(track) = track[,1]
count_tab = read.table ("C:/Users/user/Desktop/HCVandNOR/R_dada2/ASVs_counts.tsv", header=T, row.names=1, check.names=F, sep="\t" )

# Build data sheet
asv_depth = data.frame (row.names = colnames (count_tab),
                        ASV_number = colSums (count_tab != 0),  
                        depth = track["nonchim"])
colnames(asv_depth) = c("asv_number","depth")
write.csv(asv_depth,"C:/Users/user/Desktop/HCVandNOR/R_dada2/asv_depth.csv")
asv_depth$depth = as.numeric(as.character(asv_depth$depth))
asv_depth$logdepth = NA 
asv_depth$logdepth = log10(asv_depth$depth)
asv_depth$Group = c("Control","Control","Control","Control","Control","Control","Control","Control","Control","Control",
                    "Control","Control","Control","Control","Control","Control","Control","Control","Control","Control",
                    "Control","Control","Control","Control","Control","Control","Control","Control","Control","Control",
                    "Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control",
                    "SVR","SVR","SVR","BL","SVR","BL","BL","BL","BL","SVR",
                    "BL","SVR","BL","SVR","BL","SVR","BL","SVR","BL","SVR",
                    "BL","SVR","BL","SVR","BL","SVR","BL","SVR","BL","SVR",
                    "BL","SVR","BL","SVR","BL","SVR","BL","SVR","BL","SVR",
                    "BL","SVR","BL","SVR","BL","SVR","BL","SVR","BL","SVR",
                    "BL","SVR","BL","SVR","BL","SVR","BL","SVR","BL","SVR",
                    "BL","SVR","BL","SVR","BL","SVR","BL","SVR","BL","SVR",
                    "BL","SVR","BL","SVR","BL","SVR","BL","SVR","BL","SVR",
                    "SVR","BL","SVR","BL")

asv_depth$Group = factor(asv_depth$Group, levels=c("Control", "BL", "SVR"))
asv_depth$oriDep = track["input"] 

# asv_depth$oriDep = as.numeric(asv_depth$oriDep)
asv_depth$logoriDep = NA 
asv_depth$logoriDep = log10(asv_depth$oriDep)

# Draw plot
# Create a custom color scale
library(RColorBrewer)

# Asv_number/depth
myColors = brewer.pal(3,"Set1")
names(myColors) = levels(asv_depth$Group)
colScale = scale_colour_manual(name = "Group",values = myColors)
p = ggplot(asv_depth,aes(x = logdepth, y = asv_number,colour = Group)) + geom_point()
p1 = p + colScale + geom_abline()
# Legend(¡§topleft¡¨, legend = c("r = 0.32, p ¡X value = 3x10^-4","r = 0.36, p-value = 0.0"))
##+geom_smooth(method = "lm")
ggsave("C:/Users/user/Desktop/HCVandNOR/R_dada2/asv_depth_plot.png",p1)

# Calculate correlation ##here use two methods to do this
# Grab first 3 column,remove the non useful or the numeric data.
library(Hmisc)
reseq_asv_depth = read.csv("C:/Users/user/Desktop/HCVandNOR/R_dada2/reseq_asv_depth.csv")
A = as.matrix(reseq_asv_depth)[c(1:42),] # [c(43:84),],[c(85:126),]
rownames(A) = A[,1] # column 1 as rownames
A = A[,-1] # remove first column
COR = rcorr(as.matrix(A),type="pearson")
# Result:(person correlation,p-value), CON(-0.03,0.8405),BL(0.39,0.0112),SVR(0.70,0)
aA = reseq_asv_depth[c(1:42),]
rownames(aA) = aA[,1] # column 1 as rownames
aA = aA[,-1] # remove first column
aA$asv_number = as.numeric(aA$asv_number)
CORT = cor.test(aA$asv_number, aA$logdepth, method = c("pearson"))
CORT
##Result is the same as rcorr

# Boxplot of Logdepth
p = ggplot(asv_depth, aes(x = Group, y = logdepth)) +
  geom_boxplot()
ggsave("C:/Users/user/Desktop/HCVandNOR/R_dada2/depth_plot.png",p)

# Boxplot of Asv_number
p = ggplot(asv_depth, aes(x = Group, y = asv_number)) +
  geom_boxplot()
ggsave("C:/Users/user/Desktop/HCVandNOR/R_dada2/asv_plot.png",p)

# Distribution zero in each sample
prop.null = apply(count_tab, 2, function(x) 100*mean(x==0))
print(head(prop.null))
not_prop.null = 100-prop.null
not_prop.null = data.frame(not_prop.null)
not_prop.null$group = c("non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV",
                        "non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV",
                        "non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV",
                        "non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV","non-HCV",
                        "SVR","SVR","SVR","BL","SVR","BL","BL","BL","BL","SVR",
                        "BL","SVR","BL","SVR","BL","SVR","BL","SVR","BL","SVR",
                        "BL","SVR","BL","SVR","BL","SVR","BL","SVR","BL","SVR",
                        "BL","SVR","BL","SVR","BL","SVR","BL","SVR","BL","SVR",
                        "BL","SVR","BL","SVR","BL","SVR","BL","SVR","BL","SVR",
                        "BL","SVR","BL","SVR","BL","SVR","BL","SVR","BL","SVR",
                        "BL","SVR","BL","SVR","BL","SVR","BL","SVR","BL","SVR",
                        "BL","SVR","BL","SVR","BL","SVR","BL","SVR","BL","SVR",
                        "SVR","BL","SVR","BL")

not_prop.null$group <- factor(not_prop.null$group,levels=c("non-HCV","BL","SVR"))
library(RColorBrewer)
myColors = brewer.pal(3,"Set1")
names(myColors) = levels(not_prop.null$group)
colScale = scale_colour_manual(name = "group",values = myColors)
p = ggplot(not_prop.null,aes(x = group, y = not_prop.null,colour = group)) + geom_point()
p1 = p + geom_boxplot()+colScale 
ggsave("C:/Users/user/Desktop/HCVandNOR/R_dada2/not_prop.null.png",p1,dpi = 300)
dev.off()

#-------------------------------------------------------------------------------
# Setting up our working environment
#Install
## install.packages("phyloseq")
    install.packages("vegan")
## install.packages("DESeq2")
## install.packages("ggplot2")
   install.packages("dendextend")
## install.packages("tidyr")
## install.packages("viridis")
   install.packages("reshape")
#Library
  library("phyloseq")
  library("vegan")
  library("DESeq2")
  library("ggplot2")
  library("dendextend")
  library("tidyr")
  library("viridis")
  library("reshape")

# Now starting
rm(list=ls())
# Import

# Asv_table
  count_tab = read.table("C:/Users/user/Desktop/HCVandNOR/R_dada2/ASVs_counts.tsv", header=T, row.names=1, check.names=F, sep="\t")
# Taxo_table
tax_tab = as.matrix(read.table("C:/Users/user/Desktop/HCVandNOR/R_dada2/ASVs_taxonomy.tsv", header=T, row.names=1, check.names=F, sep="\t"))  # Metadata, metadata convert xlsx2tsv use:https://products.aspose.app/cells/conversion/xlsx-to-tsv
sample_info_tab = read.table("C:/Users/user/Desktop/HCVandNOR/metadata/t_metadata2.tsv", header=T, row.names=1, check.names=F, sep="\t")
# Let's add some colors to the sample info table that are specific to sample types and characteristics that we will use when plotting things
sample_info_tab = cbind(sample_info_tab,color = ("std"))
sample_info_tab$color[sample_info_tab$Cohort=="non-HCV"] = "darksalmon"
sample_info_tab$color[sample_info_tab$Cohort=="Treatment-naive"] = "chartreuse3"
sample_info_tab$color[sample_info_tab$Cohort=="SVR"] = "cornflowerblue"
# And setting the color column to be of type "character", which helps later
sample_info_tab$Cohort = as.factor ( sample_info_tab$Cohort )
sample_info_tab  #Take a peek

# Normalize each ASV length 
# Because weakly correlation between ASV length and read counts,so I dno't normalize  the ASV length.

# To get each ASV sequence length
library("Biostrings") # Read .fa file
s = readDNAStringSet("C:/Users/user/Desktop/HCVandNOR/R_dada2/ASVs.fa")
s.df = data.frame(width(s))
colnames(s.df) = c("ASV_length")

# Check correlation "pearson"
counts = data.frame(rowSums(count_tab))
cunts_length = cbind(s.df,counts)
rcorr(as.matrix(cunts_length),type="pearson") # result:-0.01
      
# Normalizing for sampling depth( VST,RLE )
# VST
count_tab_gei1 = count_tab+1 #for VST normalize because each number can't be 0
vst_deseq_matrix =  DESeqDataSetFromMatrix(count_tab_gei1, colData = sample_info_tab, design = ~ Cohort ) 
deseq_counts_vst = varianceStabilizingTransformation(vst_deseq_matrix,blind=FALSE)
# NOTE: If you get this error here with your dataset: "Error in
# estimateSizeFactorsForMatrix (counts(object), locfunc =locfunc, : every
# gene contains at least one zero, cannot compute log geometric means", that
# can be because the count table is sparse with many zeroes, which is common
# with marker-gene surveys. In that case you'd need to use a specific
# function first that is equipped to deal with that. You could run:
# deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
# now followed by the transformation function:
# deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
## For the error, I plus "one" to every counts
  
# And here is pulling out our transformed table
vst_trans_count_tab = assay(deseq_counts_vst)
    
#Convert negtive to zero
vst_trans_count_tab[vst_trans_count_tab<0]=0
    

# RLE
rle_deseq_matrix =  DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~ Cohort ) 
dds = estimateSizeFactors(rle_deseq_matrix, type = ("poscounts")) #one sample exist not only one zero
rle_trans_counts_tab = counts(dds,normalized=T) 

#--------------------------------------------------------------------------------  
# Use phyloseq
# Making our phyloseq object with transformed table
vst_count_phy = otu_table(vst_trans_count_tab, taxa_are_rows = TRUE)
vst2_count_phy = otu_table(df2, taxa_are_rows = T) #df2 is minus ASV of ASV1sample 
rle_count_phy = otu_table(rle_trans_counts_tab, taxa_are_rows = T)
sample_info_tab_phy = sample_data(sample_info_tab)
count_phy = otu_table(count_tab, taxa_are_rows = T)
count2_phy = otu_table(df22, taxa_are_rows = T)
  
# TAXO transformed table
tax_phy = tax_table(tax_tab)
# Learn from official phyloseq
physeqno = phyloseq(count_phy, tax_phy)
physeqno2 = phyloseq(count2_phy, tax_phy)
vst_physeq = phyloseq(vst_count_phy, tax_phy)
vst2_physeq = phyloseq(vst2_count_phy, tax_phy)
rle_physeq = phyloseq(rle_count_phy, tax_phy)
##plot_bar(vst_physeq, fill = "Family")
  
# Built tree   
library("ape")
vst_random_tree = rtree(ntaxa(vst_physeq), rooted = TRUE, tip.label = taxa_names(vst_physeq))
vst2_random_tree = rtree(ntaxa(vst2_physeq), rooted = TRUE, tip.label = taxa_names(vst2_physeq))
rle_random_tree = rtree(ntaxa(rle_physeq), rooted = TRUE, tip.label = taxa_names(rle_physeq))
random_treeno = rtree(ntaxa(physeqno), rooted = TRUE, tip.label = taxa_names(physeqno))
random_treeno2 = rtree(ntaxa(physeqno2), rooted = TRUE, tip.label = taxa_names(physeqno2))
##plot(vst_random_tree) #fail
  
# Merge phyloseq
m.vst_physeq = merge_phyloseq(vst_physeq, sample_info_tab_phy, vst_random_tree)
m.vst2_physeq = merge_phyloseq(vst2_physeq, sample_info_tab_phy, vst2_random_tree)
m.rle_physeq = merge_phyloseq(rle_physeq, sample_info_tab_phy, rle_random_tree)
physeqno = merge_phyloseq(physeqno, sample_info_tab_phy, random_treeno)
physeqno2 = merge_phyloseq(physeqno2, sample_info_tab_phy, random_treeno2)

# Generating and visualizing the PCoA with phyloseq
dist1 = c("unifrac")
dist2 = c("bray") #have negative number can't use bray
m.vst2_pcoa = ordinate(m.vst2_physeq, method="PCoA", distance= dist1)
m.vst_pcoa = ordinate(m.vst_physeq, method="PCoA", distance= dist1)
m.rle_pcoa = ordinate(m.rle_physeq, method="PCoA", distance= dist1)
m.rle2_pcoa = ordinate(m.rle_physeq, method="PCoA", distance= dist1)
pcoa = ordinate(physeqno, method="PCoA", distance= dist1)
pcoa2 = ordinate(physeqno2, method="PCoA", distance= dist1)

p = plot_ordination(physeqno2, pcoa2, color="Cohort") + 
  geom_point(size = 1) +
  theme(aspect.ratio = 1) 
  #geom_text(aes(label = rownames(sample_info_tab), hjust = 0.4, vjust = 0.3),check_overlap = TRUE)
  ## Three red dots mix with HCV is "CO-3289","CO-2602","CO-3457" 
  ggsave("C:/Users/user/Desktop/HCVandNOR/R_dada2/m.nonormalize_mASV1sample_unifrac_plot.png",p)
dev.off()

#-------------------------------------------------------------------------------  
# Alpha diversity ## use original count table

# Import
library(vegan)
ori_count_tab = read.table("C:/Users/user/Desktop/HCVandNOR/R_dada2/ASVs_counts.tsv", header=T, row.names=1, check.names=F, sep="\t")
ori_count_tab = t(ori_count_tab) #transpose
# Plot of relation 
S <- specnumber(ori_count_tab)
raremax <- min(rowSums(ori_count_tab))
Srare <- rarefy(ori_count_tab, raremax)
par(mfrow = c(1,2))
  
p =plot(S, Srare, xlab = "Sample size", ylab = "ASV")
   abline(0, 1)
#rarefaction curve
pdf('C:/Users/user/Desktop/HCVandNOR/R_dada2/rarecurve.pdf',width = 4.5, height = 4)
rarecurve(ori_count_tab, step = 10000, col = sample_info_tab$color,lwd=2,ylab="ASVs", xlab= "Depth", label =F, main = "Rarefaction curves")
legend("topright", legend=c("non-HCV","Treatment-naive","SVR"), fill=c("#e9967a","#66cd00","#6495ed"),cex = 0.7)  #color code
dev.off()

# First we need to create a phyloseq object using our un-transformed count table
  ori_count_tab_phy = otu_table(ori_count_tab, taxa_are_rows=T)
  tax_tab_phy = tax_table(tax_tab)
  ASV_physeq = phyloseq(ori_count_tab_phy, tax_tab_phy, sample_info_tab_phy)

# And now we can call the plot_richness() function on our phyloseq object
  sample_data(ASV_physeq)$Cohort <- factor(sample_data(ASV_physeq)$Cohort,levels=c("non-HCV","Treatment-naive","SVR"))
  p = plot_richness(ASV_physeq, "Cohort", measures=c("Simpson", "Shannon", "InvSimpson")) 
  p = p + geom_boxplot(data = p$data, aes(x = Cohort, y = value, color = NULL), 
                         alpha = 0.1)
  ggsave("C:/Users/user/Desktop/HCVandNOR/R_dada2/alphadiversity.png",p)

#--------------------------------------------------------------------------
# Boxplot beta diversity :relative abundance

#vst 
relate_vst_physeq = transform_sample_counts(m.vst_physeq, function(x) x / sum(x) * 100) 
head(otu_table(relate_vst_physeq)[,1:6]) #check

a = otu_table(relate_vst_physeq)
colSums(a[,1]) #check
disrel_unifra <- phyloseq::distance(relate_vst_physeq, method = "UniFrac")
disrel_unifra <- as.matrix(disrel_unifra)
head(disrel_unifra)[,1:6]
sub_dist <- list()
Cohort_all <- sample_data(relate_vst_physeq)$Cohort
for (Cohort in levels(Cohort_all)) { 
  row_group <- which(Cohort_all == Cohort)
  sample_group <- sample_names(relate_vst_physeq)[row_group]
  sub_dist[[Cohort]] <- disrel_unifra[ sample_group, sample_group]
  sub_dist[[Cohort]][!lower.tri(sub_dist[[Cohort]])] <- NA
}
  
unifracgroups<- melt(sub_dist)
df.unifrac <- unifracgroups[complete.cases(unifracgroups), ]
df.unifrac$L1 <- factor(df.unifrac$L1, levels=names(sub_dist))
  
head(df.unifrac)

# Place by specific order of plot
df.unifrac$L1 <- factor(sample_data(df.unifrac)$L1,levels=c("non-HCV","Treatment-naive","SVR"))
ggsave("C:/Users/user/Desktop/HCVandNOR/R_dada2/vst_beta_diversity_boxplot.png")
p= ggplot(df.unifrac, aes(x=L1, y=value, colour=L1)) +
   geom_jitter(width= 0.3) + 
   geom_boxplot(alpha=0.4,width = 0.6) +  
   theme(legend.position="none") +
   ylab("Unifrac diversity") +
   theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30,hjust=1,vjust=1,size=14), axis.text.y=element_text(size=14))+
   scale_y_continuous(limits = c(0.3,1))

# Rle        
relate_rle_physeq = transform_sample_counts(m.rle_physeq, function(x) x / sum(x) * 100) 
head(otu_table(relate_rle_physeq)[,1:6]) #check
  
a = otu_table(relate_rle_physeq)
colSums(a[,1]) #check
disrel_unifra <- phyloseq::distance(relate_rle_physeq, method = "UniFrac")
disrel_unifra <- as.matrix(disrel_unifra)
head(disrel_unifra)[,1:6]
sub_dist <- list()
Cohort_all <- sample_data(relate_rle_physeq)$Cohort
for (Cohort in levels(Cohort_all)) { 
  row_group <- which(Cohort_all == Cohort)
  sample_group <- sample_names(relate_rle_physeq)[row_group]
  sub_dist[[Cohort]] <- disrel_unifra[ sample_group, sample_group]
  sub_dist[[Cohort]][!lower.tri(sub_dist[[Cohort]])] <- NA
}
  
unifracgroups<- melt(sub_dist)
df.unifrac <- unifracgroups[complete.cases(unifracgroups), ]
df.unifrac$L1 <- factor(df.unifrac$L1, levels=names(sub_dist))
head(df.unifrac)

# Place by specific order of plot
df.unifrac$L1 <- factor(sample_data(df.unifrac)$L1,levels=c("non-HCV","Treatment-naive","SVR"))
ggsave("C:/Users/user/Desktop/HCVandNOR/R_dada2/rle_beta_diversity_boxplot.png")
p= ggplot(df.unifrac, aes(x=L1, y=value, colour=L1)) +
   geom_jitter(width= 0.3) + 
   geom_boxplot(alpha=0.4,width = 0.6) +  
   theme(legend.position="none") +
   ylab("Unifrac diversity") +
   theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30,hjust=1,vjust=1,size=14), axis.text.y=element_text(size=14))+  
   scale_y_continuous(limits = c(0.3,1)) #set range of y axis show  

#------------------------------------------------------------------------------  
#Content:1.Prevalence 2.
#prevalence hist
#column= sample
#row = ASV
#if is not zero -> 1
#prvalence = (rowSum)/126
#draw his(accumulate)
#non_chim = readRDS("C:/Users/user/Desktop/HCVandNOR/R_dada2/merge_seqtab.nochim.rds")
df = count_tab  
df[df!=0] = 1
prevalence = data.frame(prevalence = rowSums(df)/126,row.names = row.names(count_tab))
prevalence$prevalence = as.numeric(prevalence$prevalence)
png("C:/Users/user/Desktop/HCVandNOR/R_dada2/hist.prevalence.png")
hist(prevalence$prevalence, main = "Prevalence of ASV",xlab="prevalence",ylab = "frequency"  )
dev.off()
  
s_prevalence = sort.DataFrame(prevalence)
png("C:/Users/user/Desktop/HCVandNOR/R_dada2/plot.cumulative.prevalence.png")
plot.ecdf(prevalence$prevalence, main = "Prevalence of ASV",xlab="prevalence",ylab = "cumulative frequency" )
dev.off()
  
#sort my data
#0.007936508(1-21793)ASV1sample
#0.015873016(21794-24668)ASV2sample
#0.02380952(24669-26226)ASV3sample
or.prevalence = prevalence[order(prevalence$prevalence),,drop=FALSE]
  
#Get name of ASV only show in one sample
ASV1sample = subset(or.prevalence, prevalence<0.01,select = ("prevalence"))
#save name as a list
name_asv1sample = list(row.names(ASV1sample))
#
df2 = vst_trans_count_tab[-which(rownames(vst_trans_count_tab) %in% rownames(ASV1sample)), ]
df22 = count_tab[-which(rownames(count_tab) %in% rownames(ASV1sample)), ]
#USE df2
df2[df2!=0] = 1
prevalence2 = data.frame(prevalence = rowSums(df)/126,row.names = row.names(count_tab))
prevalence2$prevalence = as.numeric(prevalence2$prevalence)
png("C:/Users/user/Desktop/HCVandNOR/R_dada2/hist.prevalence2.png")
hist(prevalence2$prevalence, main = "Prevalence2 of ASV",xlab="prevalence",ylab = "frequency"  )
dev.off()

s2_prevalence = sort.DataFrame(prevalence2)
png("C:/Users/user/Desktop/HCVandNOR/R_dada2/plot.cumulative.prevalence2.png")
plot.ecdf(prevalence2$prevalence, main = "Prevalence2 of ASV",xlab="prevalence",ylab = "cumulative frequency" )
dev.off()

# Taxo remove low prevalence
m.tax_tab = tax_tab[-which(rownames(tax_tab)%in% rownames(ASV1sample)),]
dim(m.tax_tab) #9668
# Check vst table negtive value
# a = data.frame(vst_trans_count_tab[vst_trans_count_tab<0]) ,result:3872709 obs
# sum(count_tab==0) , >3872709
# All the negtive value in "vst_trans_count_tab" is from 0
  
# Check taxo_tab if there have phylum is NA, other is not NA.
# Checktaxo = tax_tab[which(is.na(tax_tab[,"Phylum"]) & !is.na(tax_tab[,"Class"]) & !is.na(tax_tab[,"Kingdom"])),]
# Result :every levels I have tried. As phylum is NA , other levels will be NA.
  
#Check (1)ASV low prevalence (2)taxonomy->phylum no identified are same groups.
library("dplyr")
l_freqiency = cbind(tax_tab,prevalence)
l_freqiency = l_freqiency %>% arrange(prevalence)
naphylum = data.frame(which (is.na(l_freqiency$Phylum)))
  
# Taxonomy
p = plot_bar(relate_vst_physeq, fill = "Genus") +
    geom_bar(aes(color = Species, fill = Species),stat = "identity",  position = "stack") +
    labs(x = "", y = "Relative Abundance\n") +
    facet_wrap(~ Cohort, scales = "free") +
    theme(panel.background = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text = element_text(size = 20))
  
  ggsave('C:/Users/user/Desktop/HCVandNOR/R_dada2/taxocomposition_genus3.png',width = 120, height = 100,limitsize = FALSE, dpi = 300,p)
  dev.off()
  
# Taxo_genus 
## Remove low prevalence ASV (only one sample exist)fail
## Could try top20 or top50 ASV
df = vst_trans_count_tab[-which(rownames(vst_trans_count_tab) %in% rownames(ASV1sample)), ]
vst2_count_phy = otu_table(df2, taxa_are_rows = T) #df2 is minus ASV of ASV1sample 
sample_info_tab_phy = sample_data(sample_info_tab)
m.tax_phy = tax_table(m.tax_tab)
vst2_physeq = phyloseq(vst2_count_phy, m.tax_phy)

# Built tree
library("ape")
vst2_random_tree = rtree(ntaxa(vst2_physeq), rooted = TRUE, tip.label = taxa_names(vst2_physeq))
vst_ASV1sample_physeq = merge_phyloseq(vst2_physeq, sample_info_tab_phy, vst2_random_tree)
relate_vst_ASV1sample_physeq = transform_sample_counts(vst_ASV1sample_physeq, function(x) x / sum(x) * 100) 
  
p = plot_bar(relate_vst_ASV1sample_physeq, fill = "Genus") +
    geom_bar(aes(color = Species, fill = Species),stat = "identity",  position = "stack") +
    labs(x = "", y = "Relative Abundance\n") +
    facet_wrap(~ Cohort, scales = "free") +
    theme(panel.background = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text = element_text(size = 20))
  
  ggsave('C:/Users/user/Desktop/HCVandNOR/R_dada2/taxocomposition_genus3.png',width = 20, height = 15,p)
  dev.off()