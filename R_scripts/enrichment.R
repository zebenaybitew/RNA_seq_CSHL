###########################################################
# Script: get_ensembl_ids.R
# Purpose: Retrieve Ensembl Gene IDs + run enrichment analysis
###########################################################

library(biomaRt)
library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)


setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/CSHL")

###########################################################
### 1) Load files
###########################################################

df <- fread("/Users/cr23646/Desktop/background genes.txt")

prot_291 <- fread(
  "/Users/cr23646/Desktop/CRC/chap_1/11_final_291_proteins.txt",
  sep = "\t", quote = FALSE, header = FALSE
)

prot_subset <- prot_291 %>%
  select(UniProt = V2, OlinkID = V3, ProteinName = V4) %>%
  mutate(across(everything(), ~ gsub('"', "", .x)))

###########################################################
### 2) Query Ensembl
###########################################################

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

results_uniprot <- getBM(
  attributes = c("uniprotswissprot", "ensembl_gene_id", "hgnc_symbol"),
  filters = "uniprotswissprot",
  values = df$UniProt,
  mart = ensembl
)

merged_data <- merge(
  df, results_uniprot,
  by.x = "UniProt", by.y = "uniprotswissprot",
  all.x = TRUE
)

fwrite(merged_data, "ensembl_ids_mapped.txt", sep = "\t", quote = FALSE)

df3 <- fread("ensembl_ids_mapped.txt")

merged_data <- prot_subset %>% left_join(df3, by = c("UniProt", "OlinkID"))

###########################################################
### 3) Background + significant lists
###########################################################

Olink_proteins <- df3 %>% select(ensembl_gene_id)
sig_proteins  <- merged_data %>% select(ensembl_gene_id)


fwrite(Olink_proteins, "Olink_proteins.txt", col.names = FALSE, quote = FALSE)
fwrite(sig_proteins,  "sig_proteins2.txt",  col.names = FALSE, quote = FALSE)


###########################################################
### 4) Convert to ENTREZ
###########################################################

sig_ids <- bitr(sig_proteins$ensembl_gene_id, fromType="ENSEMBL",
                toType="ENTREZID", OrgDb="org.Hs.eg.db")

bg_ids  <- bitr(Olink_proteins$ensembl_gene_id, fromType="ENSEMBL",
                toType="ENTREZID", OrgDb="org.Hs.eg.db")

###########################################################
### 5) Helper functions (NO repetition)
###########################################################

create_outdir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

safe_plot <- function(p, file, out) {
  if (is.null(p)) return(NULL)
  png(file.path(out, paste0(file, ".png")), 1200, 1000, res = 150)
  print(p)
  dev.off()
}

safe_table <- function(x, name, out) {
  df <- try(as.data.frame(x), silent = TRUE)
  if (!inherits(df, "try-error") && nrow(df) > 0) {
    write.csv(df, file.path(out, paste0(name, ".csv")), row.names = FALSE)
  }
}

###########################################################
### 6) General Enrichment Function
###########################################################

run_enrichment <- function(entrez, bg, outdir, label, show_n = 20) {
  
  create_outdir(outdir)
  
  go_onts <- c("BP", "MF", "CC")
  
  go_res <- lapply(go_onts, function(ont) {
    enrichGO(gene = entrez,
             OrgDb = org.Hs.eg.db,
             keyType = "ENTREZID",
             ont = ont,
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.2,
             universe = bg,
             readable = TRUE)
  })
  names(go_res) <- paste0("GO_", go_onts)
  
  ekegg <- enrichKEGG(gene = entrez, organism = "hsa",
                      keyType = "ncbi-geneid",
                      pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                      universe = bg)
  
  er <- enrichPathway(gene = entrez, universe = bg,
                      pAdjustMethod = "BH", qvalueCutoff = 0.05)
  
  # Save tables
  lapply(names(go_res), function(n)
    safe_table(go_res[[n]], paste0(label, "_", n), outdir))
  
  safe_table(ekegg, paste0(label, "_KEGG"), outdir)
  safe_table(er,    paste0(label, "_Reactome"), outdir)
  
  # Plots
  plotlist <- list(
    GO_BP = barplot(go_res$GO_BP, showCategory = show_n) + theme_minimal(),
    GO_MF = barplot(go_res$GO_MF, showCategory = show_n) + theme_minimal(),
    GO_CC = barplot(go_res$GO_CC, showCategory = show_n) + theme_minimal(),
    KEGG  = dotplot(ekegg, showCategory = show_n) + theme_minimal(),
    Reactome = dotplot(er, showCategory = show_n) + theme_minimal()
  )
  
  lapply(names(plotlist), function(n)
    safe_plot(plotlist[[n]], paste0(label, "_", n), outdir))
  
  return(list(GO = go_res, KEGG = ekegg, Reactome = er))
}

###########################################################
### 7) ENRICHMENT RUNS (ONE BLOCK ONLY)
###########################################################

results_all <- run_enrichment(sig_ids$ENTREZID, bg_ids$ENTREZID,
                              outdir = "enrichment_results/all",
                              label = "ALL")

###########################################################
### 8) Direction-specific enrichment
###########################################################

protein_df <- fread("/Users/cr23646/Desktop/CRC/chap_1/mapping/gene_mapping/mapping_GRCh38.p14_olink.txt")
protein_df <- merge(protein_df, prot_subset, by.x = "Target", by.y = "ProteinName")

protein_df <- protein_df[!duplicated(protein_df$OlinkID)]

protein_up   <- protein_df[regulation == "Upregulated"]
protein_down <- protein_df[regulation == "Downregulated"]

# Run UP
results_up <- run_enrichment(
  protein_up$entrezgene_id,
  bg_ids$ENTREZID,
  outdir = "enrichment_results/up",
  label = "UP"
)

# Run DOWN
results_down <- run_enrichment(
  protein_down$entrezgene_id,
  bg_ids$ENTREZID,
  outdir = "enrichment_results/down",
  label = "DOWN"
)



#

proteins_up <- protein_df_up[, .(entrezgene_id)]
proteins_down <- protein_df_down[, .(entrezgene_id)]


merged_up<-merge(proteins_up, bg_ids, by.x="entrezgene_id", by.y="ENTREZID", all=FALSE)

proteins_up[, entrezgene_id := as.character(entrezgene_id)]

merged_up <- merge(
  proteins_up, bg_ids,
  by.x = "entrezgene_id", by.y = "ENTREZID",
  all = FALSE
)



proteins_down[, entrezgene_id := as.character(entrezgene_id)]

merged_down <- merge(
  proteins_down, bg_ids,
  by.x = "entrezgene_id", by.y = "ENTREZID",
  all = FALSE
)


merged_down_unique <- unique(merged_down, by = "entrezgene_id")
merged_down_unique<-as.data.frame(merged_down_unique) 

merged_up_unique<-unique(merged_up, by = "entrezgene_id")
merged_up_unique<-as.data.frame(merged_up_unique) 

bg_ids_unique<-unique(bg_ids, by = "ENTREZID")
bg_ids_unique<-as.data.frame(bg_ids_unique) 

bg_ids_unique2 <- as.data.frame(bg_ids_unique[, 1])

sig_ids_unique<-unique(sig_ids, by = "ENTREZID")
sig_ids_unique<-as.data.frame(sig_ids_unique) 


# save unique proteins lists 

fwrite(bg_ids_unique, "all_proteins.txt", quote = FALSE, sep = "\t")
fwrite(bg_ids_unique2, "all_proteins1.txt", quote = FALSE, sep = "\t", col.names = FALSE)
fwrite(sig_ids_unique, "sig_proteins.txt", quote = FALSE, sep = "\t")
fwrite(merged_up_unique, "upregulated_proteins.txt", quote = FALSE, sep = "\t")
fwrite(merged_down_unique, "downregulated_proteins.txt", quote = FALSE, sep = "\t")





Olink_proteins <- df3[, .(ensembl_gene_id)]
sig_proteins  <- merged_data[, .(ensembl_gene_id)]

protein_up2   <- prot_sub[regulation == "Upregulated"]
protein_down2 <- prot_sub[regulation == "Downregulated"]
merged_up<-merge(protein_up2, df3, by.x="Assay", by.y="gene_name", all=FALSE)
merged_down<-merge(protein_down2, df3, by.x="Assay", by.y="gene_name", all=FALSE)

# Keep only the ensembl_gene_id column
merged_up <- merged_up[, .(ensembl_gene_id)]
merged_down <- merged_down[, .(ensembl_gene_id)]

# save the results 
fwrite(merged_up, "upregulated_proteins.txt", quote = FALSE, sep = "\t")
fwrite(merged_down, "downregulated_proteins.txt", quote = FALSE, sep = "\t")


