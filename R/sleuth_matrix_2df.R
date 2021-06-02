suppressMessages({
  library('cowplot')
  library('sleuth')
  library(rhdf5)
})


sample_id <- c("NAS_T-R-1","NAS_T-R-4","NAS_T-R-5","NAS_T-R-6","NAS_T-R-7","NAS_T-R-8","NAS_T-R-9","NAS_T-R-10",
               "NAS_T-R-11","NAS_T-R-12","NAS_T-R-13","NAS_T-R-14","NAS_T-R-15","NAS_T-R-16","NAS_T-R-17",
               "NAS_T-R-19","NAS_T-R-20","NAS_T-R-25","NAS_T-R-29","NAS_T-R-30","NAS_T-R-31","NAS_T-R-32",
               "NAS_T-R-33","NAS_T-R-34","NAS_T-R-35","NAS_T-R-36","NAS_T-R-37","NAS_T-R-40","NAS_T-R-41",
               "NAS_T-R-42")

filename <- dir("KallistoOutput/")
str <- strsplit(filename, "_")

kal_dirs <- file.path("G:/Mi unidad/LNC_puri/KallistoOutput",sample_id, "abundance.h5")
kal_dirs

condition <- c(1,0,1,1,1,1,0,0,0,1,0,1,1,1,1,0,0,0,1,1,0,0,0,1,1,0,0,1,1,1)



metadata <- data.frame(sample_id, condition, kal_dirs)
names(metadata) <- c("sample","condition","path")
metadata

load("ttg_homo_sapiens_ona.RData")

for(t in 1:3)
{
  metadata[,t] <- as.character(metadata[,t])
}


# Analysis

so <- sleuth_prep(metadata, target_mapping = ttg,
                  aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE)

# obs_norm_filt = Filtered and normalized values; obs_norm = Normalized values; obs_raw = Raw values

summary(so)
names(so)
head(so$obs_raw)
dim(so$obs_raw)
head(so$obs_norm)
dim(so$obs_norm)
summary(so$obs_norm$est_counts[so$obs_norm$sample == "NAS_T-R-1"])

head(so$filter_df)
dim(so$filter_df)

head(so$obs_norm_filt)
dim(so$obs_norm_filt)
summary(so$obs_norm_filt$est_counts[so$obs_norm_filt$sample == "NAS_T-R-1"])

so$sample_to_covariates
so$target_mapping
so$gene_column

# ("obs_norm" or "obs_raw")
# ("tpm" or "est_counts" (for transcript-level analyses) or "scaled_reads_per_base" (for gene-level analyses))

sleuth_matrix_norm_tpm <- sleuth_to_matrix(so, 'obs_norm', 'est_counts')
summary(sleuth_matrix_norm_tpm)
dim(sleuth_matrix_norm_tpm)
par(mfrow = c(2,1))
boxplot(log(sleuth_matrix_norm_tpm+0.5), main = "obs_norm, est_counts")
# save(so, file = "20210223_So_kallisto.RData")
summary(so$obs_norm_filt$est_counts)
length(unique((so$obs_norm_filt$target_id)))
hist(so$obs_norm_filt$est_counts, breaks = 100000, xlim = c(0,1000))
boxplot(log(so$obs_norm_filt$est_counts+0.5) ~ so$obs_norm_filt$sample, main = "obs_norm_filt, estimcounts")

par(mfrow = c(1,1))
boxplot(log(so$obs_raw$est_counts+0.5) ~ so$obs_raw$sample, main = "obs_raw, estimcounts")

# Paciente 41

# PCA ---------------------------------------------------------------------

library(reshape2)
mat <- reshape2::acast(so$obs_norm_filt[,1:3], target_id ~ sample)



library(factoextra)

res.pca <- prcomp(t(mat), scale = TRUE)
fviz_eig(res.pca)

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

kol <- ifelse(so$sample_to_covariates$condition == 0, "black","red")

fviz_pca_ind(res.pca,
             axes = c(1, 2),
             col.ind = kol, # Color by the quality of representation
             repel = TRUE,     # Avoid text overlapping
             title = "NASIR - Progesión", addEllipses = T
)
