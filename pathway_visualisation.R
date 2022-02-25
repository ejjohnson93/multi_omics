#### Map fold change onto specific pathways ####

# The variable 'data' could be the output of a Limma or EdgeR differential expression analysis containing differentially expressed genes and LogFC values
# Pathway enrichment analysis should have already been carried out so pathways of interest are already known 

# Load libaries

library(SBGNview)
library(KEGG.db)
library(annotate)
library(org.Hs.eg.db)
library(pathview)


#### Pathview ####

fold_changes <- dplyr::select(data, "logFC")

# Annotate with Entrez IDs instead (needed for pathview)
cols <- c("ENTREZID", "SYMBOL")
annotation <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(fold_changes), columns=cols, keytype="SYMBOL")
fold_changes <- cbind(fold_changes, entrez = annotation[,2])

# Remove repeat IDS
fold_changes <- fold_changes[!duplicated(fold_changes$entrez), ]
fold_changes <- fold_changes[!is.na(fold_changes$entrez), ]
rownames(fold_changes) <- fold_changes$entrez

fold_changes <- dplyr::select(fold_changes, "logFC")

# Example calcium signalling pathway
pv_out <- pathview(gene.data = fold_changes, pathway.id = "hsa04020",
                   species = "hsa", out.suffix = "fleischer", kegg.native = T, limit = list(gene = 0.01, cpd = 1))



#### SBGNview ####

# Load all pathway data in SBGN database
data("sbgn.xmls")

#If needed convert gene IDs
#Example to convert entrez IDs into pathCommons IDs as SBGNview has functionality to do so
gene_data <- changeDataId(data.input.id = fold_changes,
                          input.type = "entrez",
                          output.type = "pathwayCommons",
                          mol.type = "gene",
                          sum.method = "sum")

# Pathways of interest can be queried using key words
# In this case I just searched the name of the pathway but multiple search terms could be provided

findPathways(c("Crosslinking of collagen fibrils"))

# Create visualisation for collagen cross-linking pathway
SBGNview_obj <- SBGNview(gene.data = gene_data,
                         input.sbgn = "R-HSA-2243919",
                         max.gene.value = 0.1,
                         min.gene.value = -0.1,
                         mid.gene.value = 0,
                         output.file = "collagen_cross-linking",
                         gene.id.type = "pathwayCommons",
                         output.formats =  c("png", "pdf", "ps"))
SBGNview_obj


