#### Load Libraries ####
#BiocManager::install(organism, character.only = TRUE)
#install.packages("ggridges")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(tidyverse)
library(ggridges)

#### Set Organism ####
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

#### Load & Tidy Data ####
# reading in log fold change data
df = read_csv("Fred_Hutch_R/data/24_naive_lfc_data.csv")
colnames(df)[1] = "gene"
colnames(df)[2] = "lfc"

# we want the log2 fold change
original_gene_list <- df$lfc

# name the vector
names(original_gene_list) <- df$gene

# omit any NA values
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

#### gseGO ####
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

#### Visualizations ####
# Dotplot
require(DOSE)
dotplot(gse, showCategory=10, split=".sign", font.size = 8) + facet_grid(.~.sign)

# Enrichment Map
gse_emap_data <- pairwise_termsim(gse)
emapplot(gse_emap_data, showCategory = 10)

# Category Netplot
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, categorySize="pvalue", color.params = list(foldChange = gene_list), showCategory = 3)

# Ridge Plot
ridgeplot(gse) + labs(x = "enrichment distribution")
