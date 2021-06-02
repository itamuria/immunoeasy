suppressMessages({
  library('cowplot')
  library('sleuth')
  library(rhdf5)

  library(officer)
  library(flextable)
  library(tidyverse)
  library(tidyquant)
  library(timetk)
})

setwd("G:/Mi unidad/KarmeleValencia/DSTYK/Kallisto_musmus/")

# Define variables
group <- "3LL"
go_file <- paste0("G:/Mi unidad/KarmeleValencia/DSTYK/results/20201014_GO_",group,".xlsx")
so_file <- paste0("G:/Mi unidad/KarmeleValencia/DSTYK/rdata/20201014_So_mmu_",group,".RData")
genes_file <- paste0("G:/Mi unidad/KarmeleValencia/DSTYK/results/20201014_Genes_",group,".xlsx")
##



#  Analysis ---------------------------------------------------------------



sample_id <- c("3LL_empty_R1_S7",
               "3LL_empty_R2_S8",
               "3LL_empty_R3_S9",

               "3LL_sh3_R1_S10",
               "3LL_sh3_R2_S11",
               "3LL_sh3_R3_S12")



filename <- dir()
str <- strsplit(filename, "_")

kal_dirs <- file.path(sample_id, "abundance.h5")
kal_dirs

condition <- c(rep("Empty",3), rep("Sh3",3))

metadata <- data.frame(sample_id, condition, kal_dirs)
names(metadata) <- c("sample","condition","path")
metadata

load("ttg.RData")

for(t in 1:3)
{
  metadata[,t] <- as.character(metadata[,t])
}


# Analysis

so <- sleuth_prep(metadata, target_mapping = ttg,
                  aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE)



so <- sleuth_fit(so, ~ 1, 'reduced')
so <- sleuth_fit(so, ~ condition, 'full')
so <- sleuth_lrt(so, 'reduced', 'full')

save(so, file = so_file)
# load(so_file)

pca_plot <- plot_pca(so, color_by = 'condition')
pca_plot2 <- plot_pca(so, color_by = 'condition', text_labels = TRUE)


sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_table_gene <- dplyr::filter(sleuth_table_gene, qval <= 0.05)
dim(sleuth_table_gene)
# head(sleuth_table_gene, 20)

# transcript
sleuth_table_tx <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = FALSE)
sleuth_table_tx <- dplyr::filter(sleuth_table_tx, qval <= 0.05)
dim(sleuth_table_tx)
# head(sleuth_table_tx, 20)

# Wald test for a sleuth model
oe <- sleuth_wt(so, which_beta = 'conditionSh3')
sleuth_results_oe <- sleuth_results(oe,
                                    test = 'conditionSh3',
                                    show_all = TRUE)

sleuth_results_oe_disagr <- sleuth_results(oe, test = 'conditionSh3',
                                           test_type = "wt", which_model = "full",
                                           pval_aggregate = FALSE)

library(dplyr)
sleuth_results_oe_disagr <- sleuth_results_oe_disagr[order(sleuth_results_oe_disagr$pval),]
sig_transcripts <- sleuth_results_oe_disagr %>% filter(qval < 0.05)
dim(sig_transcripts)
# plot_transcript_heatmap(oe,
# transcripts = sig_transcripts$target_id[1:50])


library(EnhancedVolcano)

volc <- EnhancedVolcano(sleuth_results_oe_disagr,
                        lab = sleuth_results_oe_disagr$ext_gene,
                        x = 'b',
                        y = 'pval',
                        # xlim = c(-8, 8),
                        title = paste0(group, ' empty vs. sh3'),
                        pCutoff = 10e-40,
                        FCcutoff = 1,
                        pointSize = 3.0,
                        labSize = 3.0)


# Plotting

dstk_plot <- plot_bootstrap(oe,
                            target_id = "ENSMUST00000045110.13",
                            units = "est_counts",
                            color_by = "condition")


plot_bootstrap(oe,
               target_id = "ENSMUST00000069949.12",
               units = "est_counts",
               color_by = "condition")

# How many with different qval

library(dplyr)

temp <- sleuth_results_oe_disagr %>% filter(qval < 0.05)
dim(temp)
length(unique(temp$ext_gene))


# Funcionalities ----------------------------------------------------------

openxlsx::write.xlsx(temp, file = genes_file)

# Load libraries
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Mm.eg.db)

## Run GO enrichment analysis
ego <- enrichGO(gene = temp$ens_gene,
                universe = sleuth_results_oe_disagr$ens_gene,
                keyType = "ENSEMBL",
                OrgDb = org.Mm.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

## Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

# write.csv(cluster_summary, "results/clusterProfiler_Mov10oe.csv")
openxlsx::write.xlsx(cluster_summary, file = go_file)

## Dotplot
dotplot <- dotplot(ego, showCategory=50)

## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
emaplot <- emapplot(ego, showCategory = 50)


## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
OE_foldchanges <- temp$b

names(OE_foldchanges) <- temp$ext_gene

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cneplot <- cnetplot(ego,
                    categorySize="pvalue",
                    showCategory = 5,
                    foldChange=OE_foldchanges,
                    vertex.label.font=6)


## Run GO enrichment analysis
ego100 <- enrichGO(gene = temp$ens_gene[1:100],
                   universe = sleuth_results_oe_disagr$ens_gene,
                   keyType = "ENSEMBL",
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)

## Output results from GO analysis to a table
cluster_summary100 <- data.frame(ego100)

# write.csv(cluster_summary, "results/clusterProfiler_Mov10oe.csv")
go_file100 <- paste0("G:/Mi unidad/KarmeleValencia/DSTYK/results/20201014_GO_",group,"100.xlsx")
openxlsx::write.xlsx(cluster_summary100, file = go_file100)


## If some of the high fold changes are getting drowned out due to a large range, you could set a maximum fold change value
# OE_foldchanges <- ifelse(OE_foldchanges > 2, 2, OE_foldchanges)
# OE_foldchanges <- ifelse(OE_foldchanges < -2, -2, OE_foldchanges)
#
# cnetplot(ego,
#          categorySize="pvalue",
#          showCategory = 5,
#          foldChange=OE_foldchanges,
#          vertex.label.font=6)

# Power point -------------------------------------------------------------


doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(doc, value = volc, location = ph_location_fullsize())
doc <- ph_with(doc, value = "Volcano plot", location = ph_location_type(type = "title"))

doc <- add_slide(doc)
doc <- ph_with(doc, value = pca_plot, location = ph_location_fullsize())

doc <- add_slide(doc)
doc <- ph_with(doc, value = pca_plot2, location = ph_location_fullsize())

doc <- add_slide(doc)
doc <- ph_with(doc, value = dstk_plot, location = ph_location_fullsize())

doc <- add_slide(doc)
doc <- ph_with(doc, value = dotplot, location = ph_location_fullsize())

doc <- add_slide(doc)
doc <- ph_with(doc, value = emaplot, location = ph_location_fullsize())

doc <- add_slide(doc)
doc <- ph_with(doc, value = cneplot, location = ph_location_fullsize())


# doc <- ph_with(doc, value = volc, location = ph_location_left())
# doc <- ph_with(doc, value = stock_plot, location = ph_location_right())

print(doc, target = paste0("G:/Mi unidad/KarmeleValencia/DSTYK/results/dstyk_results_",group,".pptx"))

# print results

warning(paste0("Results of ",group))
warning(dim(temp))
warning(length(unique(temp$ext_gene)))
warning(dim(sleuth_table_tx)[1])
