# QUANTILE normalization --------------------------------------------------

setwd("G:/Mi unidad/NASIR/SignaturesSandra/")

load("G:/Mi unidad/NASIR/NASIR_counts/NASIR_subread_allCounts_annotated_fpmk2.RData")

head(annotLookup2)
ann <- annotLookup2
pat <- openxlsx::read.xlsx("G:/Mi unidad/NASIR/MAIN_FOLDER_NASIR/MAIN_TABLE_v002.xlsx",2)
names(pat)
pat2 <- pat[,c(1,10)]
pat3 <- pat2[!is.na(pat2$TUMOR_RNA),]

names(ann)
ann2 <- ann[,c(2,64:93)]
head(ann2)

colnames(ann2) <- c("symbol", pat3$TUMOR_RNA)
ann3 <- ann2[ann2$symbol != "",]
ann4 <- data.frame(t(ann3[,-1]))
names(ann4) <- ann3$symbol

# save(ann4, file = "ann4.RData")



# Quitar los que no tengan datos ni varianza ------------------------------

tf_sum0 <- apply(ann4, 2, function(x) sum(x, na.rm = T)) == 0
ann5 <- ann4[,!tf_sum0]

# No hay var == 0
# tf_var0 <- apply(ann5, 2, function(x) var(x, na.rm = T)) == 0
# ann6 <- ann5[,!tf_var0]

# Quitar con algún NA

tf_sum_NA <- apply(ann5, 2, function(x) sum(is.na(x), na.rm = T)) > 0
ann6 <- ann5[,!tf_sum_NA]

# quantile ----------------------------------------------------------------


# library(raster)
# x <- 1:100

quantile(ann5$TSPAN6, seq(0, 10, by = 0.1))

library(gtools)
decil <- quantcut(ann4$TSPAN6, seq(0,1,by=0.1) )
levels(decil)

h <- 0
ann7 <- sapply(ann6, function(x)
{
  h <<- h + 1
  print(paste0(h, " - ", names(ann5)[h]))
  decil <- quantcut(x, seq(0,1,by=0.1) )
  as.numeric(decil)
})

ann7b <- data.frame(Pat = rownames(ann6), ann7)


# Log10 -------------------------------------------------------------------

ann8_log <- log10(ann7b[,-1])
ann8_log <- data.frame(Pat = rownames(ann6), ann8_log)

save(ann8_log, file = "ann8_log.RData")

# Groups mean -----------------------------------------------------------

#• Load groups of genes

gengr <- openxlsx::read.xlsx("G:/Mi unidad/NASIR/SignaturesSandra/Gene-signature.xlsx",3)

# Check if we have all the interesting genes

gengr$Genes[!gengr$Genes %in% colnames(ann7)]

un_gr <- unique(gengr$Group)

un_gen <- unique(gengr$Genes)
un_gen <- un_gen[which(un_gen != "NGK7")]

colnames(ann8_log) <- gsub("-",".",colnames(ann8_log) )
un_gen <- gsub("-",".",un_gen)

grep("HLA",colnames(ann8_log), value = T)

library(tidyverse)
ann9_fil <- tibble(data.frame(Pat = rownames(ann6),ann8_log[,colnames(ann8_log) %in% un_gen]))

library(dplyr)

names(ann9_fil)
grep("HLA",names(ann9_fil), value = T)

ann9_fil %>%
  mutate(Inflamatory = CD274+CD8A+LAG3+STAT1,
         Cytolytic32 = GZMA+PRF1,
         Gajewski33 = CCL2+CCL3+CCL4+CD8A+CXCL10+CXCL9+GZMK+HLA.DMA+HLA.DMB+HLA.DOA+HLA.DOB+ICOS+IRF1,
         InterferonGamaSig = CXCL10+CXCL9+HLA.DRA+IDO1+IFNG+STAT1,
         AntigenPresenting = CMKLR1+HLA.DQA1+HLA.DRB1+PSMB10,
         InterferonGamaBiology = CCL5+CD27+CXCL9+CXCR6+IDO1+STAT1,
         TcellExhaustion = CD276+CD8A+LAG3+PDCD1LG2+TIGIT,
         T_nK_sig = HLA.E,
         RibasGeneInterferon = CCR5+CXCL10+CXCL11+CXCL9+GZMA+HLA.DRA+IDO1+IFNG+PRF1+STAT1,
         Inflamatory_mean = (CD274+CD8A+LAG3+STAT1)/4,
         Cytolytic32_mean = (GZMA+PRF1)/2,
         Gajewski33_mean = (CCL2+CCL3+CCL4+CD8A+CXCL10+CXCL9+GZMK+HLA.DMA+HLA.DMB+HLA.DOA+HLA.DOB+ICOS+IRF1)/13,
         InterferonGamaSig_mean = (CXCL10+CXCL9+HLA.DRA+IDO1+IFNG+STAT1)/6,
         AntigenPresenting_mean = (CMKLR1+HLA.DQA1+HLA.DRB1+PSMB10)/4,
         InterferonGamaBiology_mean = (CCL5+CD27+CXCL9+CXCR6+IDO1+STAT1)/6,
         TcellExhaustion_mean = (CD276+CD8A+LAG3+PDCD1LG2+TIGIT)/5,
         T_nK_sig_mean = (HLA.E)/1,
         RibasGeneInterferon_mean = (CCR5+CXCL10+CXCL11+CXCL9+GZMA+HLA.DRA+IDO1+IFNG+PRF1+STAT1)/10) -> ann10_fin

ann9_fil %>% select(CD274,CD8A,LAG3,STAT1) %>%
  transmute(Inflamatory = CD274 + CD8A + LAG3 + STAT1)
