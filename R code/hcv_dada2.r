#@For 16s HCV data,DADA2 v1.8
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
  path <- "/home/blair/HCV_CON/HCV" 
  ##list.files(path)
  # Forward and reverse fastq filenames have format: SAMPLENAME_R*.clean.fastq.gz
  fnFs <- sort(list.files(path, pattern="_R1.clean.fastq.gz", full.names = TRUE))
  fnRs <- sort(list.files(path, pattern="_R2.clean.fastq.gz", full.names = TRUE))
  # Extract sample names,ex:DAA-E041-BL_R2.clean.fastq.gz
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  saveRDS(sample.names,"/home/blair/HCV_CON/R_dada2/samplenames/hcv_sams.rds")
  ## Inspect read quality profiles
  # forward
  f = plotQualityProfile(fnFs[1:2])
  ggsave("/home/blair/HCV_CON/R_dada2/QC_beforetrunc/hcvquality_fwd.png",f)
  dev.off()
  # reverse
  r = plotQualityProfile(fnRs[1:2])
  ggsave("/home/blair/HCV_CON/R_dada2/QC_beforetrunc/hcvquality_rvs.png",r)
  dev.off()
  
  
# Filter and trim
# place filtered files in filtered/ subdirectory
  filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(255,240),
                       trimLeft=c(17,21),trimRight=c(5,30),maxN=0, maxEE=c(3,5), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=TRUE) 
  write.csv(out,file="/home/blair/HCV_CON/R_dada2/truncresult/hcvtruncresult.csv")
  ###forward
  a = plotQualityProfile(filtFs[1:2])
  ggsave("/home/blair/HCV_CON/R_dada2/QC_aftertrunc/hcvtrunc_quality_fwd.png",a)
  dev.off()
  ###reverse
  b = plotQualityProfile(filtRs[1:2])
  ggsave("/home/blair/HCV_CON/R_dada2/QC_aftertrunc/hcvtrunc_quality_rvs.png",b)
  dev.off()
  out
  
# Learn the Error Rates
  errF <- learnErrors(filtFs, multithread=TRUE)
  c = plotErrors(errF, nominalQ=TRUE)
  ggsave("/home/blair/HCV_CON/R_dada2/errorrates_plot/hcverrorrates_fwd.png",c)
  dev.off()
  saveRDS(errF,file = "/home/blair/HCV_CON/R_dada2/errRDS/hcverrF.rds")
  
  errR <- learnErrors(filtRs, multithread=TRUE)
  d = plotErrors(errR, nominalQ=TRUE)
  ggsave("/home/blair/HCV_CON/R_dada2/errorrates_plot/hcverrorrates_rwd.png",d)
  dev.off()
  saveRDS(errR,file = "/home/blair/HCV_CON/R_dada2/errRDS/hcverrR.rds")

# Dereplicate
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  saveRDS(derepFs,"/home/blair/HCV_CON/R_dada2/derepRDS/derepFs.rds")
  
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  saveRDS(derepRs,"/home/blair/HCV_CON/R_dada2/derepRDS/derepRs.rds")
# Name the derep-class objects by the sample names
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names

  
# Sample Inference
  dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
  saveRDS(dadaFs,"/home/blair/HCV_CON/R_dada2/dadaRDS/dadaFs.rds")
  dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
  saveRDS(dadaRs,"/home/blair/HCV_CON/R_dada2/dadaRDS/dadaRs.rds")
  ##dadaFs[[1]]
  
  
# Merge paired reads , minoverlap default is 12bp.my merge result is not good so I loose the rule.
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 10, maxMismatch = 0, verbose=TRUE)
  # Inspect the merger data.frame from the first sample
  head(mergers[[1]])
  saveRDS(mergers,file="/home/blair/HCV_CON/R_dada2/mergersRDS/hcv_mergers.rds")


#Construct sequence table
  seqtab <- makeSequenceTable(mergers)
  ## The sequences being tabled vary in length.
  dim(seqtab)
  ## Inspect distribution of sequence lengths
  # I don't save this result for non-hav group,because I have viewed it already when processing.  
  a = table(nchar(getSequences(seqtab)))
  write.table(a, file="/home/blair/HCV_CON/hcv_lengthdisttribution.csv")

    
# Remove chimeras, and construct sequence table
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", minFoldParentOverAbundance = 8,multithread=TRUE, verbose=TRUE)
  dim(seqtab.nochim)
  sum(seqtab.nochim)/sum(seqtab)
  saveRDS(seqtab.nochim, "/home/blair/HCV_CON/R_dada2/seqtab.nochimRDS/hcvseqtab.nochim.rds")

  
# Track reads through the pipeline
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
  # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names
  head(track)
  write.csv(track,file="/home/blair/HCV_CON/R_dada2/track/hcvtrack.csv")

#DADA2 finished!  
  