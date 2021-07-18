#@For 16s non-HCV data,DADA2 v1.8
#@
#@The data pre-procedure,and generating sequence table by DADA2
#@
#@author:tzu-yu Chu
#@data:06.29.2021

#@-----------------------------------------------


library("dada2")
library("ggplot2")
# Import data
# file contain HCV fastq.gz
path <- "/home/blair/HCV_CON/Control" 
##list.files(path)
# For change non-hcv sample names
names = 2426
# Forward and reverse fastq filenames have format: SAMPLENAME_R*.clean.fastq.gz
fnFs <- sort(list.files(path, pattern=paste0(names,"_R1.clean.fastq.gz"), full.names = TRUE))
fnRs <- sort(list.files(path, pattern=paste0(names,"_R2.clean.fastq.gz"), full.names = TRUE))
# Extract sample names,ex:DAA-E041-BL_R2.clean.fastq.gz
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
saveRDS(sample.names,paste0("/home/blair/HCV_CON/R_dada2/samplenames/",names,"_sams.rds"))
## Inspect read quality profiles
# forward
f = plotQualityProfile(fnFs[1:2])
ggsave(paste0("/home/blair/HCV_CON/R_dada2/QC_beforetrunc/",names,"quality_fwd.png"),f)
dev.off()
# reverse
r = plotQualityProfile(fnRs[1:2])
ggsave(paste0("/home/blair/HCV_CON/R_dada2/QC_beforetrunc/",names,"quality_rvs.png"),r)
dev.off()


# Filter and trim
# place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(245,255),
                     trimLeft=c(21,21),trimRight=c(25,25),maxN=0, maxEE=c(3,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
out
write.csv(out,file=paste0("/home/blair/HCV_CON/R_dada2/truncresult/",names,"truncresult.csv"))
###forward
f = plotQualityProfile(filtFs[1:2])
ggsave(paste0("/home/blair/HCV_CON/R_dada2/QC_aftertrunc/",names,"trunc_quality_fwd.png"),f)
dev.off()
###reverse
r = plotQualityProfile(filtRs[1:2])
ggsave(paste0("/home/blair/HCV_CON/R_dada2/QC_aftertrunc/",names,"trunc_quality_rvs.png"),r)
dev.off()


# Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
f = plotErrors(errF, nominalQ=TRUE)
ggsave(paste0("/home/blair/HCV_CON/R_dada2/errorrates_plot/",names,"errorrates_fwd.png"),f)
dev.off()
saveRDS(errF,file = paste0("/home/blair/HCV_CON/R_dada2/errRDS/",names,"errF.rds"))

errR <- learnErrors(filtRs, multithread=TRUE)
r = plotErrors(errR, nominalQ=TRUE)
ggsave(paste0("/home/blair/HCV_CON/R_dada2/errorrates_plot/",names,"errorrates_rwd.png"),r)
dev.off()
saveRDS(errR,file = paste0("/home/blair/HCV_CON/R_dada2/errRDS/",names,"errR.rds"))


# Dereplicate
derepFs <- derepFastq(filtFs, verbose=TRUE)
saveRDS(derepFs,paste0("/home/blair/HCV_CON/R_dada2/derepRDS/",names,"derepFs.rds"))

derepRs <- derepFastq(filtRs, verbose=TRUE)
saveRDS(derepRs,paste0("/home/blair/HCV_CON/R_dada2/derepRDS/",names,"derepRs.rds"))


# Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
saveRDS(dadaFs,paste0("/home/blair/HCV_CON/R_dada2/dadaRDS/",names,"dadaFs.rds"))
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
saveRDS(dadaRs,paste0("/home/blair/HCV_CON/R_dada2/dadaRDS/",names,"dadaRs.rds"))
##dadaFs[[1]]


# Merge paired reads , minoverlap default is 12bp.my merge result is not good so I loose the rule.
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 12, maxMismatch = 0, verbose=TRUE)
  # Inspect the merger data.frame from the first sample
  head(mergers[[1]])
  saveRDS(mergers,file=paste0("/home/blair/HCV_CON/R_dada2/mergersRDS/",names,"_mergers.rds"))


#Construct sequence table
  seqtab <- makeSequenceTable(mergers)
  ## The sequences being tabled vary in length.
  dim(seqtab)
  # Inspect distribution of sequence lengths
  table(nchar(getSequences(seqtab)))


# Remove chimeras, and construct sequence table
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", minFoldParentOverAbundance = 8,multithread=TRUE, verbose=TRUE)
  dim(seqtab.nochim)
  sum(seqtab.nochim)/sum(seqtab)
  rownames(seqtab.nochim) = sample.names
  saveRDS(seqtab.nochim, paste0("/home/blair/HCV_CON/R_dada2/seqtab.nochimRDS/",names,"seqtab.nochim.rds"))


# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, getN(dadaFs), getN(dadaRs), getN(mergers), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track,file=paste0("/home/blair/HCV_CON/R_dada2/track/",names,"track.csv"))

#DADA2 finished!




