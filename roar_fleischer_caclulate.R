######################
## Install and load packages ----

#install.packages("BiocManager")
#BiocManager::install("tidyverse")
#BiocManager::install("roar")
#BiocManager::install("biomaRt")
#BiocManager::install("tximport")


library(tidyverse)
library(roar)
library(biomaRt)
library(tximport)


## Read in samples ----

# Lists of paths to BAM files 
bamTreatment <- list.files(path = "/mnt/lustre/users/ejohn16/RNAseq-files/datasets/human/ageing_datasets/skin_fleischer_PRJNA454681/hisat2",
                           pattern = ".bam",
                           full.names = TRUE)

bamControl <- list.files(path = "/mnt/lustre/users/ejohn16/RNAseq-files/datasets/human/ageing_datasets/skin_fleischer_PRJNA454681/hisat2",
                         pattern = ".bam",
                         full.names = TRUE)


# Extract individual file names 
sample_names <- list.files(path = "/mnt/lustre/users/ejohn16/RNAseq-files/datasets/human/ageing_datasets/skin_fleischer_PRJNA454681/hisat2",
                         pattern = ".bam",
                         full.names = FALSE)

sample_names <- sample_names %>% gsub(pattern = "_hisat2.bam", replacement = "", x = .) 

names(bamTreatment) <- paste0(sample_names, "_treatment")
names(bamControl) <- paste0(sample_names, "_control")


# GTF file path ----

gtf <- "/mnt/data1/users/ejohn16/RNAseq/genomes/human/hg19/GTF/hs_apasdb.gtf"




## Create RoarDataset objects ----
# Separately for each BAM file


allRes <- list()
allfpkm <- list()

for(i in 1:length(bamTreatment)){
    print(i)
    # Create roar dataset object
    rds <- RoarDatasetFromFiles(bamTreatment[i], 
                                          bamControl[i], 
                                          gtf)
    
    # Calculate ratios 
    rds <- countPrePost(rds, FALSE)
    rds <- computeRoars(rds)
    
    # Get ratios and fpkm 
    res <- fpkmResults(rds)
    allRes[[i]] <- res %>% dplyr::select("mM_treatment")
    allfpkm[[i]] <- res %>% dplyr::select("treatmentValue")
   
}

# Bind together ratios for all samples, remove columns for treatment
allRes <- dplyr::bind_cols(allRes)
allfpkm <- dplyr::bind_cols(allfpkm)



## Process raw data ----
# name columns 

metadata <- read.csv("/mnt/lustre/users/ejohn16/RNAseq-files/datasets/human/ageing_datasets/skin_fleischer_PRJNA454681/roar/fleischer_metadata.csv", header = TRUE, stringsAsFactors = F)

print(sample_names == metadata$Run)

colnames(allRes) <- t(paste0(metadata$Run, "_ratio"))
colnames(allfpkm) <- t(paste0(metadata$Run, "_fpkm"))


# remove NAs
# 'Inf' is used for NA values in roar so these need to be reassigned as NA first
allRes[allRes == Inf] <- NA
allRes <- allRes[rowSums(is.na(allRes)) != ncol(allRes),]




## Annotate and save fpkm data (to be used later) ----
geneID_fpkm <- rownames(allfpkm) %>% gsub(pattern = "_.*",
                                    replacement = "",
                                    x = .)

allfpkm <- allfpkm %>% mutate(entrezgene = as.numeric(geneID_fpkm))

write_csv(allfpkm, "/mnt/lustre/users/ejohn16/RNAseq-files/datasets/human/ageing_datasets/skin_fleischer_PRJNA454681/roar/fleischer_FPKM_table.csv")



## Get annotation ----

geneID <- rownames(allRes) %>% gsub(pattern = "_.*",
                                    replacement = "",
                                    x = .)

allRes <- allRes %>% mutate(entrezgene = as.numeric(geneID))


# get gene names from gene id, hg19

mart <- useEnsembl(biomart = "ensembl",
                      dataset = "hsapiens_gene_ensembl",
                      host="http://apr2019.archive.ensembl.org")


annotation <- getBM(attributes = c("entrezgene","external_gene_name"), 
                    filters = "entrezgene", 
                    values = geneID, 
                    ensembl,
                    mart = mart)

results_anno <- left_join(allRes, annotation)




## Get TPMs ----

files <- list.files(path = "/mnt/lustre/users/ejohn16/RNAseq-files/datasets/human/ageing_datasets/skin_fleischer_PRJNA454681/salmon_quant/",
                    pattern = "quant.sf",
                    recursive = TRUE,
                    full.names = TRUE)

names(files) <- t(paste0(metadata$Run, "_TPM"))

# create tx2gene table
transcript_ID <- read.table(files[1], sep = "\t", header = TRUE)[,1] %>%
                    as.character(.) %>%
                    gsub(pattern = "\\.\\d*", 
                         replacement = "", 
                         x = .)

annotation2 <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", 
                                    "external_gene_name", "description",
                                    "gene_biotype"), 
                    filters = "ensembl_transcript_id", 
                    values = transcript_ID, 
                    mart = mart)

transcript_ID_original <- read.table(files[1], 
                                     sep = "\t", 
                                     header = TRUE)[,1] %>%
                            as.character(.) %>% 
                            data.frame(transcript_versionID = ., transcript_ID)

tx2gene_key <- merge(transcript_ID_original, annotation2[,c(1:2)], 
                     by.x = "transcript_ID",
                     by.y = "ensembl_transcript_id") %>% .[,(2:3)]


txi.salmon <- tximport(files, 
                       type = "salmon", 
                       tx2gene = tx2gene_key, 
                       countsFromAbundance = "lengthScaledTPM")


# table of TPM values for all samples
TPM_table <- txi.salmon$abundance %>%
    data.frame(.) %>%
    mutate(ensembl_gene_id = row.names(.)) %>%
    left_join(., annotation2[,c(2:5)],
              by = c("ensembl_gene_id")) %>%
    unique(.)



## Process tables and save data ----

# check sum of TPMs in each sample, should be 1 million
sum(TPM_table[,1])
            
# combined table
results_anno_TPM <- left_join(results_anno, TPM_table) 

#sum(results_anno_TPM$PRJNA223350_3_TPM, na.rm = TRUE)

# save results - tables
write.table(results_anno, file = "/mnt/lustre/users/ejohn16/RNAseq-files/datasets/human/ageing_datasets/skin_fleischer_PRJNA454681/roar/fleischer_results_anno.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)


write.table(results_anno_TPM, file = "/mnt/lustre/users/ejohn16/RNAseq-files/datasets/human/ageing_datasets/skin_fleischer_PRJNA454681/roar/fleischer_sep_roar_results.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

write.table(TPM_table, file = "/mnt/lustre/users/ejohn16/RNAseq-files/datasets/human/ageing_datasets/skin_fleischer_PRJNA454681/roar/fleischer_TPM_table.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)


# save results - csv files
write_csv(results_anno_TPM, "/mnt/lustre/users/ejohn16/RNAseq-files/datasets/human/ageing_datasets/skin_fleischer_PRJNA454681/roar/fleischer_sep_roar_results.csv")

write_csv(TPM_table, "/mnt/lustre/users/ejohn16/RNAseq-files/datasets/human/ageing_datasets/skin_fleischer_PRJNA454681/roar/fleischer_TPM_table.csv")
#save.image("~/Desktop/Emily_APA/Datasets/PRJNA223350/ROAR/Working_ROAR_Workspace.RData")