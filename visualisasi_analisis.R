##Analisis awal (konfirmasi keberadaan gen signifikan)
#set directory
hasil_kegg <- read.csv("./proteins_KEGG_enrichment.csv")
gene_map <- read.csv("./genemapped.csv")

#menggabungkan hasil gene_map dengan kegg berdasarkan protein.id
merged_data <- merge(hasil_kegg, gene_map, by = "Protein.Id", all.x = TRUE)

#konfirmasi keberhasilan hasil merge
head(merged_data)


#Loading package yang digunakan 
library(tidyverse)
"dplyr" %in% .packages()

### konfirmasi protein signifikan
#Daftar gen yang berkaitan dengan penyakit NASH
daftar_gen_nash <- c("Srebf1", "Ppara", "Pparg", "Tnfa", "Tgfbr2", "Mapk8", 
                    "Nfkb1", "Pik3r1", "Akt1", "Mtor", "Rxra", "Fasn", "Prkaa1")

#filtering protein
protein_nash <- merged_data %>%
  filter(Gene.x %in% daftar_gen_nash) %>%
  distinct(Protein.Id, .keep_all = TRUE) # Menghilangkan redundansi protein

#Hasil keberadaan protein terkait
print(protein_nash[, c("Protein.Name.x", "Gene.x", "Pathway_Name")])


### RANK protein
#Loading package yang digunakan 
library(STRINGdb)
library(igraph)
library(tidyverse)

#filtering spesifik menggunakan package dplyr
protein_nash <- merged_data %>%
  dplyr::filter(Gene.x %in% daftar_gen_nash) %>%
  dplyr::distinct(Protein.Id, .keep_all = TRUE)

# Mengkalkulasi interaksi protein
string_db <- STRINGdb$new(version = "12.0", species = 10090, score_threshold = 400)
mapped_proteins <- string_db$map(protein_nash, "Gene.x", removeUnmappedRows = TRUE)
interaksi_ppi <- string_db$get_interactions(mapped_proteins$STRING_id)

#Perhitungan PPI
hitung_ppi <- graph_from_data_frame(interaksi_ppi, directed = FALSE)

#Pemeringkatan interaksi dan keterkaitan protein
centrality_scores <- tibble(
  STRING_id = V(hitung_ppi)$name,
  Degree = degree(hitung_ppi, v = V(hitung_ppi), mode = "all"),
  Betweenness = betweenness(hitung_ppi, v = V(hitung_ppi), directed = FALSE)
) %>%
  dplyr::left_join(mapped_proteins %>% dplyr::select(STRING_id, Gene.x), by = "STRING_id") %>%
  dplyr::select(Gene = Gene.x, Degree, Betweenness) %>%
  dplyr::distinct() %>%
  arrange(desc(Betweenness), desc(Degree))

print("10 Teratas protein NASH dengan interaksi tertinggi")
print(head(centrality_scores, 10))


### KEGG PATHWAY
#loading package yang akan digunakan 
library(clusterProfiler)
library(org.Mm.eg.db) # database mouse
library(ggplot2)


#Pencarian protein unik pada data NASH yang sudah dibuat sebelumnya
gen_terkait_nash <- merged_data %>%
  filter(Gene.x %in% daftar_gen_nash) %>%
  pull(Gene.x) %>% # Dijadikan vektor
  unique()         # Memastikan protein unik

#Analisis KEGG
entrez_ids <- bitr(gen_terkait_nash, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

#enrichment analysis
hasil_enrichment <- enrichKEGG(
  gene = entrez_ids$ENTREZID,
  organism = 'mmu', # 'mmu' adalah kode mouse
  pvalueCutoff = 0.05
)

#visualisasi barplot
barplot(hasil_enrichment, showCategory = 15) +
  labs(title = "Pathway KEGG yang berkaitan dengan protein NASH")


##PPI (Protein-Protein Interaction)

#Loading package yang digunakan 
  library(STRINGdb)
  
#Pendefinisian stringdb
string_db <- STRINGdb$new(
  version = "12.0",
  species = 10090,  # CRITICAL FIX: 10090 is for Mus musculus (mouse)
  score_threshold = 400, # A medium confidence score
  network_type = "full"
)

#filtering protein
protein_nash <- merged_data %>%
  filter(Gene.x %in% daftar_gen_nash) %>%
  distinct(Protein.Id, .keep_all = TRUE)

#mapping terhadap database string
mapped_proteins <- string_db$map(protein_nash, "Gene.x", removeUnmappedRows = TRUE)

#pembuatan plot ppi
string_db$plot_network(mapped_proteins$STRING_id)
