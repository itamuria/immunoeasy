load("G:/Mi unidad/PabloSarobe/Immunodeconv/20210121_CIBERSORT_Results.RData")

rownames(results) <- c(colnames(tt)[2:14],pat3$ID)

write.table(results[-c(1:13),], file = "20210201_CIBERSORT_results.txt", sep = "\t", row.names = TRUE, quote = FALSE)
res2 <- data.frame(results[-c(1:13),])
res2$Pat <- rownames(res2)
openxlsx::write.xlsx(res2, file = "20210201_CIBERSORT_results.xlsx")

setwd("G:/Mi unidad/PabloSarobe/Immunodeconv/")

library(Rserve)
library('e1071')
library('parallel')
library('colorRamps')

Rserve(args="--no-save")

source("G:/Mi unidad/PabloSarobe/Immunodeconv/CIBERSORT.R")

lm22 <- read.table("LM22.txt", sep = "\t")

results <- CIBERSORT("LM22.txt", "20210121_symbol_fpkm_nasir.txt", perm=100, QN=TRUE)
