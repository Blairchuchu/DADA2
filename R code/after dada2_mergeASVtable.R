#@For 16s HCV, non-HCV seperate table
#@
#@The merge and build ASV table
#@
#@author:tzu-yu Chu
#@data:06.31.2021

#@-----------------------------------------------
library("dada2")
library("ggplot2")

k2p= readRDS("C:/Users/user/Desktop/HCVandNOR/mergetable_taxo.k2p_silva_nr99_v138.rds")

# Rename sample names (on asv table and I also rename metadata in excel.Let then have the same ID names for each samples.)
asv_table= read.table("C:/Users/user/Desktop/HCVandNOR/ASVs_counts.tsv")
write.table(colnames(asv_table),"C:/Users/user/Desktop/HCVandNOR/columnnames.csv")
colnames(asv_table)=
  c("CO-2087",
    "CO-2091",
    "CO-2116",
    "CO-2348",
    "CO-2355",
    "CO-2370",
    "CO-2371",
    "CO-2388",
    "CO-2395",
    "CO-2398",
    "CO-2399",
    "CO-2414",
    "CO-2423",
    "CO-2426",
    "CO-2429",
    "CO-2438",
    "CO-2495",
    "CO-2519",
    "CO-2543",
    "CO-2560",
    "CO-2602",
    "CO-2653",
    "CO-2677",
    "CO-2684",
    "CO-2699",
    "CO-2828",
    "CO-2916",
    "CO-3178",
    "CO-3195",
    "CO-3215",
    "CO-3231",
    "CO-3242",
    "CO-3255",
    "CO-3271",
    "CO-3275",
    "CO-3289",
    "CO-3457",
    "CO-3470",
    "CO-3473",
    "CO-3479",
    "CO-3493",
    "CO-3506",
    "102-DAA-0075-SVR",
    "102-DAA-0077-SVR",
    "102-DAA-0078-SVR",
    "102-DAA0070-BL",
    "102-DAA0070-SVR",
    "102-DAA0075-BL",
    "102-DAA0077-BL",
    "102-DAA0078-BL",
    "102-DAA0080-BL",
    "102-DAA0080-SVR",
    "102-DAA0081-BL",
    "102-DAA0081-SVR",
    "102-DAA0084-BL",
    "102-DAA0084-SVR",
    "102-DAA0086-BL",
    "102-DAA0086-SVR",
    "102-DAA0087-BL",
    "102-DAA0087-SVR",
    "102-DAA0088-BL",
    "102-DAA0088-SVR",
    "102-DAA0089-BL",
    "102-DAA0089-SVR",
    "102-DAA0091-BL",
    "102-DAA0091-SVR",
    "102-DAA0095-BL",
    "102-DAA0095-SVR",
    "102-DAA0096-BL",
    "102-DAA0096-SVR",
    "102-DAA0097-BL",
    "102-DAA0097-SVR",
    "102-DAA0098-BL",
    "102-DAA0098-SVR",
    "102-DAA0099-BL",
    "102-DAA0099-SVR",
    "103-DAA-E045-BL",
    "103-DAA-E045-SVR",
    "103-DAA-E047-BL",
    "103-DAA-E047-SVR",
    "103-DAA-E048-BL",
    "103-DAA-E048-SVR",
    "103-DAA-E050-BL",
    "103-DAA-E050-SVR",
    "103-DAA-E051-BL",
    "103-DAA-E051-SVR",
    "103-DAA-E053-BL",
    "103-DAA-E053-SVR",
    "103-DAA-E057-BL",
    "103-DAA-E057-SVR",
    "103-DAA-E059-BL",
    "103-DAA-E059-SVR",
    "103-DAA-E061-BL",
    "103-DAA-E061SVR",
    "103-DAA-E062-BL",
    "103-DAA-E062-SVR",
    "103-DAA-E063-BL",
    "103-DAA-E063-SVR",
    "103-DAA-E066-BL",
    "103-DAA-E066-SVR",
    "DAA-E010-BL",
    "DAA-E010-SVR",
    "DAA-E012-BL",
    "DAA-E012-SVR",
    "DAA-E014-BL",
    "DAA-E014-SVR",
    "DAA-E015-BL",
    "DAA-E015-SVR",
    "DAA-E016-BL",
    "DAA-E016-SVR",
    "DAA-E020-BL",
    "DAA-E020-SVR",
    "DAA-E021-BL",
    "DAA-E021-SVR",
    "DAA-E023-BL",
    "DAA-E023-SVR",
    "DAA-E038-BL",
    "DAA-E038-SVR",
    "DAA-E041-BL",
    "DAA-E041-SVR",
    "DAA-E043-BL",
    "DAA-E043-SVR",
    "DAA0046-SVR",
    "DAA0046-BL",
    "DAA0060-SVR",
    "DAA0060-BL")
  write.table(asv_table,"C:/Users/user/Desktop/HCVandNOR/ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)


# Merge hcv& non-hcv nochim tables (on Ubuntu)
  # List path
    path <- "/home/blair/HCV_CON/R_dada2/seqtab.nochimRDS" 
    ##list.files(path)
  # Forward and reverse fastq filenames have format: SAMPLENAME_R*.clean.fastq.gz
    seqtabpath <- list.files(path, pattern=".nochim.rds", full.names = TRUE)
    write.table(seqtabpath,"/home/blair/HCV_CON/R_dada2/seqtab.nochimRDS/pathlist_seqtab.nochim.txt")
  # Merge HCV and non-HCV tables
    st1 = readRDS("/home/blair/HCV_CON/R_dada2/seqtab.nochimRDS/2087seqtab.nochim.rds")
    st2 = readRDS("/home/blair/HCV_CON/R_dada2/seqtab.nochimRDS/2091seqtab.nochim.rds")
    mergeseqtable = mergeSequenceTables(st1,st2,st3,st4,st5,st6,st7,st8,st9,st10,
                                        st11,st12,st13,st14,st15,st16,st17,st18,st19,st20,
                                        st21,st22,st23,st24,st25,st26,st27,st28,st29,st30,
                                        st31,st32,st33,st34,st35,st36,st37,st38,st39,st40,
                                        st41,st42,st43)
    saveRDS(mergeseqtable,file ="/home/blair/HCV_CON/R_dada2/seqtab.nochimRDS/merge_seqtab.nochim.rds") 


# Taxo (on Ubuntu)
  # silva_nr99_v138 Kimdom to Genus  
  # wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
  # silva_nr99_v138 Species  
  # wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz
    mergetable = readRDS("/home/blair/HCV_CON/R_dada2/seqtab.nochimRDS/merge_seqtab.nochim.rds")
    taxa = assignTaxonomy(mergetable, "/home/blair/HCV_CON/R_dada2/SILVA_nr99_v138train2021.RDate", multithread=TRUE)
    saveRDS(taxa,file = "/home/blair/HCV_CON/R_dada2/mergetable_taxo_silva_nr99_v138.rds")
  # Taxo add species
    taxa = readRDS("/home/blair/HCV_CON/R_dada2/mergetable_taxo_silva_nr99_v138.rds")
    taxa.k2p <- addSpecies(taxa, "/home/blair/HCV_CON/R_dada2/SILVA_nr99_train_species/silva_species_assignment_v138.1.fa.gz", verbose=TRUE)
    saveRDS(taxa.k2p,file = "/home/blair/HCV_CON/R_dada2/mergetable_taxo.k2p_silva_nr99_v138.rds")

    
# Generate asv-table(counts) and taxo-table
  mergetable = readRDS("C:/Users/user/Desktop/HCVandNOR/R_dada2/merge_seqtab.nochim.rds")
  taxa.k2p = readRDS("C:/Users/user/Desktop/HCVandCON/R_dada2/mergetable_taxo.k2p_silva_nr99_v138.rds")
  # Giving our seq headers more manageable names (ASV_1, ASV_2...)
    asv_seqs <- colnames(mergetable)
    asv_headers <- vector(dim(mergetable)[2], mode="character")
  
    for (i in 1:dim(mergetable)[2]) {
      asv_headers[i] <- paste(">ASV", i, sep="_")
    }
  # Making and writing out a fasta of our final ASV seqs:
    asv_fasta <- c(rbind(asv_headers, asv_seqs))
    ##write(asv_fasta, "/home/blair/HCV_CON/R_dada2/ASVs.fa")
  # Count table:
    asv_tab <- t(mergetable)
    row.names(asv_tab) <- sub(">", "", asv_headers)
    ##write.table(asv_tab, "/home/blair/HCV_CON/R_dada2/ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)
  
  # Tax table:
    asv_tax = taxa.k2p
    rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)
    ##write.table(asv_tax, "/home/blair/HCV_CON/R_dada2/ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)
