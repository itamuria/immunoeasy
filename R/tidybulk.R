

suppressMessages({
  library('cowplot')
  library('sleuth')
  library(rhdf5)
})


sample_id <- c("NAS_T-R-1","NAS_T-R-4","NAS_T-R-5","NAS_T-R-6","NAS_T-R-7","NAS_T-R-8","NAS_T-R-9","NAS_T-R-10",
               "NAS_T-R-11","NAS_T-R-12","NAS_T-R-13","NAS_T-R-14","NAS_T-R-15","NAS_T-R-16","NAS_T-R-17",
               "NAS_T-R-19","NAS_T-R-20","NAS_T-R-25","NAS_T-R-29","NAS_T-R-30","NAS_T-R-31","NAS_T-R-32",
               "NAS_T-R-33","NAS_T-R-34","NAS_T-R-35","NAS_T-R-36","NAS_T-R-37","NAS_T-R-40","NAS_T-R-41",
               "NAS_T-R-42","Hepa_10619","Hepa_10622","Hepa_10627","Hepa_10628","Hepa_10632","Hepa_10634",
               "Hepa_10635","Hepa_10584","Hepa_10594","Hepa_10615","Hepa14","Hepa19","Hepa20","Hepa21",
               "Hepa23","Hepa26","CUN_10396","CUN_10426","CUN_10446","CUN_10454","CUN_10474","CUN_10509",
               "CUN_10525","CUN_10532","CUN_10543","CUN_10610","CUN_10629","CUN_10639","CUN_10640",
               "CUN_10642","CUN_10646","CUN_10647")

filename <- dir("KallistoOutput/")
str <- strsplit(filename, "_")

kal_dirs <- file.path("KallistoOutput",sample_id, "abundance.h5")
kal_dirs

condition <- c(rep("Nasir",30), rep("Hepa",16), rep("CUN",16))

metadata <- data.frame(sample_id, condition, kal_dirs)
names(metadata) <- c("sample","condition","path")
metadata

load("G:/Mi unidad/NASIR/ttg_homo_sapiens_ona.RData")

for(t in 1:3)
{
  metadata[,t] <- as.character(metadata[,t])
}




mer <- merge(fil, ttg, by = "target_id", all.x =  TRUE)
mer2 <- merge(mer, so$sample_to_covariates, by = "sample", all.x = T)

kountk <- as.tibble(mer2[,c("sample","ext_gene","est_counts","condition")])
names(kountk)[2:3] <- c("transcript","count")




library(tidybulk)

tt = kountk %>% tidybulk(sample, transcript, count)

tt %>%
  # call analysis functions
  get_bibliography()

tt.aggr = tt %>% aggregate_duplicates()
tt.norm = tt.aggr %>% identify_abundant(factor_of_interest = condition) %>% scale_abundance()

tt.norm %>%
  ggplot(aes(count_scaled + 1, group=sample, color=condition)) +
  geom_density() +
  scale_x_log10()

tt.norm %>%
  ggplot(aes(count_scaled + 1, group=sample, color=sample)) +
  geom_density() +
  scale_x_log10()
# +
#   my_theme

tt.norm.variable = tt.norm %>% keep_variable()

tt.norm.MDS =
  tt.norm %>%
  reduce_dimensions(method="MDS", .dims = 6)

tt.norm.MDS %>% pivot_sample()  %>% select(contains("Dim"), everything())

tt.norm.MDS %>%
  pivot_sample() %>%
  GGally::ggpairs(columns = 6:9, ggplot2::aes(colour=condition))

tt.norm.PCA =
  tt.norm %>%
  reduce_dimensions(method="PCA", .dims = 6)

tt.norm.PCA %>% pivot_sample() %>% select(contains("PC"), everything())

tt.norm.PCA %>%
  pivot_sample() %>%
  GGally::ggpairs(columns = 6:9, ggplot2::aes(colour=condition))

# tt.norm.tSNE =
#   tt.aggr %>%
#   identify_abundant() %>%
#   reduce_dimensions(
#     method = "tSNE",
#     perplexity=10,
#     pca_scale =TRUE
#   )
#
#
# tt.norm.tSNE %>%
#   pivot_sample() %>%
#   select(contains("tSNE"), everything())
#
# tt.norm.tSNE %>%
#   pivot_sample() %>%
#   ggplot(aes(x = `tSNE1`, y = `tSNE2`, color=Call)) + geom_point()
#
# tt.norm.MDS.rotated =
#   tt.norm.MDS %>%
#   rotate_dimensions(`Dim1`, `Dim2`, rotation_degrees = 45, action="get")
#
# tt.norm.MDS.rotated %>%
#   ggplot(aes(x=`Dim1`, y=`Dim2`, color=`Cell type` )) +
#   geom_point()
#
#
# tt.norm.MDS.rotated %>%
#   ggplot(aes(x=`Dim1 rotated 45`, y=`Dim2 rotated 45`, color=`Cell type` )) +
#   geom_point()

tt.de =
  tt.aggr %>%
  test_differential_abundance( ~ condition, action="get")
tt.de

# tt.de =
#   tt %>%
#   identify_abundant(factor_of_interest = condition) %>%
#   test_differential_abundance(
#     ~ 0 + condition,
#     .contrasts = c( "condition1 - condition0"),
#     action="get"
#   )
#
# tt.norm.adj =
#   tt.norm %>% adjust_abundance(   ~ factor_of_interest + batch)

tt.cibersort =
  tt.aggr %>%
  deconvolve_cellularity(action="get", cores=1)

library(tidyverse)
tt.cibersort %>%
  select(contains("cibersort:"), everything()) %>%
  gather(`Cell type inferred`, `proportion`, 1:22) %>%
  ggplot(aes(x=`Cell type inferred`, y=proportion, fill=condition)) +
  geom_boxplot() +
  facet_wrap(~condition) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), aspect.ratio=1/5)

tt.aggr %>%
  test_differential_cellularity( ~ condition )

# tt %>%
#   test_differential_cellularity(survival::Surv(time, dead) ~ .)

# tt.norm.cluster = tt.norm.MDS %>%
#   cluster_elements(method="kmeans", centers = 2, action="get" )
#
# tt.norm.cluster %>%
#   ggplot(aes(x=`Dim1`, y=`Dim2`, color=`cluster kmeans`)) +
#   geom_point()
#
# tt.norm.SNN =
#   tt.norm.tSNE %>%
#   cluster_elements(method = "SNN")
#
# tt.norm.SNN %>%
#   pivot_sample() %>%
#   select(contains("tSNE"), everything())
#
# tt.norm.SNN %>%
#   pivot_sample() %>%
#   gather(source, Call, c("cluster SNN", "Call")) %>%
#   distinct() %>%
#   ggplot(aes(x = `tSNE1`, y = `tSNE2`, color=Call)) + geom_point() + facet_grid(~source)

# Do differential transcription between clusters
# tt.norm.SNN %>%
#   mutate(factor_of_interest = `cluster SNN` == 3) %>%
#   test_differential_abundance(
#     ~ factor_of_interest,
#     action="get"
#   )
#
# tt.norm.non_redundant =
#   tt.norm.MDS %>%
#   remove_redundancy(    method = "correlation" )
#
# tt.norm.non_redundant %>%
#   pivot_sample() %>%
#   ggplot(aes(x=`Dim1`, y=`Dim2`, color=`Cell type`)) +
#   geom_point()

# counts_ensembl %>% ensembl_to_symbol(ens)
tt %>% describe_transcript() %>% select(transcript, description, everything())



# tidybulk ----------------------------------------------------------------
# https://github.com/stemangiola/tidybulk
# devtools::install_github("stemangiola/tidybulk")
library(tidybulk)

tt = counts %>% tidybulk(sample, transcript, count)

tt %>%
  # call analysis functions
  get_bibliography()

tt.aggr = tt %>% aggregate_duplicates()
tt.norm = tt.aggr %>% identify_abundant(factor_of_interest = condition) %>% scale_abundance()

tt.norm %>%
  ggplot(aes(count_scaled + 1, group=sample, color=`Cell type`)) +
  geom_density() +
  scale_x_log10()
# +
#   my_theme

tt.norm.variable = tt.norm %>% keep_variable()

tt.norm.MDS =
  tt.norm %>%
  reduce_dimensions(method="MDS", .dims = 6)

tt.norm.MDS %>% pivot_sample()  %>% select(contains("Dim"), everything())

tt.norm.MDS %>%
  pivot_sample() %>%
  GGally::ggpairs(columns = 10:15, ggplot2::aes(colour=`Cell type`))

tt.norm.PCA =
  tt.norm %>%
  reduce_dimensions(method="PCA", .dims = 6)

tt.norm.PCA %>% pivot_sample() %>% select(contains("PC"), everything())

tt.norm.PCA %>%
  pivot_sample() %>%
  GGally::ggpairs(columns = 10:12, ggplot2::aes(colour=`Cell type`))

tt.norm.tSNE =
  breast_tcga_mini %>%
  tidybulk(       sample, ens, count_scaled) %>%
  identify_abundant() %>%
  reduce_dimensions(
    method = "tSNE",
    perplexity=10,
    pca_scale =TRUE
  )


tt.norm.tSNE %>%
  pivot_sample() %>%
  select(contains("tSNE"), everything())

tt.norm.tSNE %>%
  pivot_sample() %>%
  ggplot(aes(x = `tSNE1`, y = `tSNE2`, color=Call)) + geom_point()

tt.norm.MDS.rotated =
  tt.norm.MDS %>%
  rotate_dimensions(`Dim1`, `Dim2`, rotation_degrees = 45, action="get")

tt.norm.MDS.rotated %>%
  ggplot(aes(x=`Dim1`, y=`Dim2`, color=`Cell type` )) +
  geom_point()


tt.norm.MDS.rotated %>%
  ggplot(aes(x=`Dim1 rotated 45`, y=`Dim2 rotated 45`, color=`Cell type` )) +
  geom_point()

tt.de =
  tt %>%
  test_differential_abundance( ~ condition, action="get")
tt.de

tt.de =
  tt %>%
  identify_abundant(factor_of_interest = condition) %>%
  test_differential_abundance(
    ~ 0 + condition,
    .contrasts = c( "conditionTRUE - conditionFALSE"),
    action="get"
  )

tt.norm.adj =
  tt.norm %>% adjust_abundance(   ~ factor_of_interest + batch)

tt.cibersort =
  tt %>%
  deconvolve_cellularity(action="get", cores=1)

library(tidyverse)
tt.cibersort %>%
  select(contains("cibersort:"), everything()) %>%
  gather(`Cell type inferred`, `proportion`, 1:22) %>%
  ggplot(aes(x=`Cell type inferred`, y=proportion, fill=`Cell type`)) +
  geom_boxplot() +
  facet_wrap(~`Cell type`) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), aspect.ratio=1/5)

tt %>%
  test_differential_cellularity( ~ condition )

tt %>%
  test_differential_cellularity(survival::Surv(time, dead) ~ .)

tt.norm.cluster = tt.norm.MDS %>%
  cluster_elements(method="kmeans", centers = 2, action="get" )

tt.norm.cluster %>%
  ggplot(aes(x=`Dim1`, y=`Dim2`, color=`cluster kmeans`)) +
  geom_point()

tt.norm.SNN =
  tt.norm.tSNE %>%
  cluster_elements(method = "SNN")

tt.norm.SNN %>%
  pivot_sample() %>%
  select(contains("tSNE"), everything())

tt.norm.SNN %>%
  pivot_sample() %>%
  gather(source, Call, c("cluster SNN", "Call")) %>%
  distinct() %>%
  ggplot(aes(x = `tSNE1`, y = `tSNE2`, color=Call)) + geom_point() + facet_grid(~source)

# Do differential transcription between clusters
tt.norm.SNN %>%
  mutate(factor_of_interest = `cluster SNN` == 3) %>%
  test_differential_abundance(
    ~ factor_of_interest,
    action="get"
  )

tt.norm.non_redundant =
  tt.norm.MDS %>%
  remove_redundancy(    method = "correlation" )

tt.norm.non_redundant %>%
  pivot_sample() %>%
  ggplot(aes(x=`Dim1`, y=`Dim2`, color=`Cell type`)) +
  geom_point()

# counts_ensembl %>% ensembl_to_symbol(ens)
tt %>% describe_transcript() %>% select(transcript, description, everything())
