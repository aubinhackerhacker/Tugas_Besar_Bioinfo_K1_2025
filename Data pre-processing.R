## Data Preprocessing
install.packages("arrow")
#Loading package
library(tidyverse)
library(arrow)

#Set directory
#Input file DIA-NN (.tsv)
data <- read.delim("/report.protein_description.tsv", header = TRUE, sep = "\t")

#Ekstraksi kolom dan pembersihan missing values
proteins_cleaned <- data %>%
  select(Protein.Name, Gene, Protein.Id) %>%
  na.omit()
write.csv(proteins_cleaned, "/bersih.csv", row.names = FALSE)


## Gene Mapping
install.packages("biomaRt")
#Loading Package
library(biomaRt)

#Input data yang telah di preprocessing
protein_df <- read.csv("./bersih.csv", stringsAsFactors = FALSE)

#Connect dengan data ensembl mouse 
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#konversi Protein ID menjadi Entrez ID and Ensembl ID
protein_info <- getBM(
  attributes = c("uniprotswissprot", "entrezgene_id", "ensembl_gene_id"),
  filters = "uniprotswissprot",
  values = unique(protein_df$Protein.Id),
  mart = ensembl
)

#Penggabungan data protein dengan info gen
final_data <- merge(protein_df, protein_info,
                    by.x = "Protein.Id", by.y = "uniprotswissprot",
                    all.x = TRUE) %>%
  filter(!is.na(entrezgene_id) & !is.na(ensembl_gene_id)) %>%
  rename(Entrez_ID = entrezgene_id,
         Ensembl_ID = ensembl_gene_id)
write.csv(final_data, "./genemapped.csv", row.names = FALSE)


## Pathway Enrichment Analysis 
# Loading package
library(clusterProfiler)
library(tidyr)

#Loading data gen yang telah di map sebelumnya
final_data <- read.csv("./genemapped.csv", stringsAsFactors = FALSE)

#KEGG pathway enrichment analysis
kegg_results <- enrichKEGG(
  gene = final_data$Entrez_ID,
  organism = "mmu",
  keyType = "kegg"
)

#Ekstraksi informasi pathway dan kalkulasi enrichment scores
kegg_df <- as.data.frame(kegg_results) %>%
  select(ID, Description, geneID, pvalue, p.adjust, GeneRatio) %>%
  mutate(Enrichment_Score = sapply(strsplit(GeneRatio, "/"),
                                   function(x) as.numeric(x[1]) / as.numeric(x[2]))) %>%
  separate_rows(geneID, sep = "/") %>%
  rename(Pathway_ID = ID,
         Pathway_Name = Description,
         Entrez_ID = geneID,
         P_Value = pvalue,
         Adjusted_P_Value = p.adjust)

#Penggabungan hasil enrichment dengan data protein-gen
final_results <- merge(final_data, kegg_df, by = "Entrez_ID", all.x = TRUE) %>%
  select(Protein.Id, Protein.Name, Gene, Entrez_ID, Ensembl_ID,
         Pathway_ID, Pathway_Name, Enrichment_Score) %>%
  distinct()
write.csv(final_results, "./proteins_KEGG_enrichment.csv", row.names = FALSE)

