#### Load Library ####
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggplot2)
library(tidyverse)

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

# If data is not already sorted: sort the list in decreasing order (required for clusterProfiler)
#gene_list = sort(gene_list, decreasing = TRUE)

#### Choose Organism ####
m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame

#### C7 gene set (Term to Gene) ####
#### in category: change to either: H, C1,C2, C3...C7
C7_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
head(C7_t2g)

#### GSEA ####
gse <- GSEA(gene_list, TERM2GENE = C7_t2g)

#### Visualizations ####
# Dotplot
require(DOSE)
dotplot(gse, showCategory=20, split=".sign", font.size = 8) + facet_grid(.~.sign)

# Enrichment Map
gse_emap_data <- pairwise_termsim(gse)
emapplot(gse_emap_data, showCategory = 50,min_edge = 0.2, cex_label_category = 0.5)

# Category Netplot
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, categorySize="pvalue", color.params = list(foldChange = gene_list), showCategory = 50)

# Ridge Plot
ridgeplot(gse, label_format = 40) + labs(x = "enrichment distribution") + theme(axis.text.y = element_text(size = 6))

# GSEA Plot

## get geneSetID number from Excel sheet of gene set (i.e Hallmark, C7,etc) to input for 
# index for  description and geneSetID
gseaplot(gse, by = "all", title = gse$Description[8], geneSetID = 8)
