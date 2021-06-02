
setwd("G:/Mi unidad/KarmeleValencia/DSTYK/TCGA2/")

# Jakiteko nola deitu
datatable(readr::read_csv("https://docs.google.com/spreadsheets/d/1f98kFdj9mxVDc1dv4xTZdx8iWgUiDYO-qiFJINvmTZs/export?format=csv&gid=2046985454",col_types = readr::cols()),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 40),
          rownames = FALSE)

# Estructura de los nombres del TCGA
# https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/


library(TCGAbiolinks)
library(dplyr)
library(DT)

# Obs: The data in the legacy database has been aligned to hg19

query.counts <-  GDCquery(project = "TCGA-LUAD",
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          workflow.type = "HTSeq - Counts")
GDCdownload(query.counts)
query.countsLUAD <- GDCprepare(query = query.counts, save = TRUE, save.filename = "LungLUAD.rda")

load("LungLUAD.rda")
query.countsLUAD <- data

query.CNV <-  GDCquery(project = "TCGA-LUAD",
                       data.category = "Copy Number Variation",
                       data.type = "Gene Level Copy Number Scores")
GDCdownload(query.CNV)
query.CNVLUAD <- GDCprepare(query = query.CNV, save = TRUE, save.filename = "../LungLUAD_CNV.rda")

# SNP
query.maf <- GDCquery(project = "TCGA-LUAD",
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Raw Simple Somatic Mutation",
                      access = "open",
                      workflow.type = "SomaticSniper Variant Aggregation and Masking")
GDCdownload(query.maf)
maf3 <- GDCprepare(query.maf)
# https://www.bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/mutation.html
# maf <- GDCquery_Maf("LUAD", pipelines = "somaticsniper") %>% read.maf
# # dmaf <- data.frame()
# datatable(maf[1:20,],
#           filter = 'top',
#           options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
#           rownames = FALSE)
#
# library(maftools)
# library(dplyr)
# datatable(getSampleSummary(maf),
#           filter = 'top',
#           options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
#           rownames = FALSE)
# plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
# oncoplot(maf = maf, top = 10, removeNonMutated = TRUE)
# titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
# #plot titv summary
# plotTiTv(res = titv)

# Frogk -------------------------------------------------------------------
library(SummarizedExperiment)

# dd <- datatable(as.data.frame(colData(query.countsLUAD)),
#                 options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
#                 rownames = FALSE)

dd <- as.data.frame(colData(query.countsLUAD))
dd2 <- as.data.frame(assay(query.countsLUAD))

table(dd$definition)
summary(dd$intermediate_dimension)
table(dd$is_ffpe)
table(dd$tissue_type)
table(dd$ajcc_pathologic_stage)
table(dd$tumor_stage)
table(dd$tissue_or_organ_of_origin)
summary(dd$days_to_last_follow_up)
table(!is.na(dd$days_to_last_follow_up))
summary(dd$age_at_diagnosis)
table(dd$primary_diagnosis)
table(dd$prior_malignancy)
table(dd$ajcc_pathologic_t)
table(dd$ajcc_pathologic_n)
table(dd$ajcc_pathologic_m)
summary(dd$cigarettes_per_day)
summary(dd$years_smoked)
summary(dd$pack_years_smoked)
table(dd$ethnicity)
table(dd$gender)
table(dd$race)
table(dd$vital_status)

dd_resumen <- dd[,c("definition","intermediate_dimension","is_ffpe","tissue_type","ajcc_pathologic_stage",
                    "tumor_stage","tissue_or_organ_of_origin",
                    "age_at_diagnosis","primary_diagnosis","prior_malignancy",
                    "ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m",
                    "cigarettes_per_day","years_smoked","pack_years_smoked","ethnicity",
                    "gender","race",
                    "days_to_last_follow_up","vital_status")]

openxlsx::write.xlsx(dd_resumen, file = "20210128_ResumenVariables.xlsx")

dataPrep_Target <- TCGAanalyze_Preprocessing(object = query.countsLUAD,
                                             cor.cut = 0.6,
                                             datatype = "HTSeq - Counts")

# https://bioconductor.org/packages/release/bioc/vignettes/EDASeq/inst/doc/EDASeq.html

dataNorm_Target <- TCGAanalyze_Normalization(tabDF = dataPrep_Target,
                                             geneInfo = geneInfoHT,
                                             method = "gcContent")

boxplot(dataPrep_Target, outline = FALSE)

boxplot(dataNorm_Target, outline = FALSE)

# https://seqqc.wordpress.com/2020/02/17/removing-low-count-genes-for-rna-seq-downstream-analysis/https://seqqc.wordpress.com/2020/02/17/removing-low-count-genes-for-rna-seq-downstream-analysis/

dataNorm <- TCGAanalyze_Normalization(tabDF = dataNorm_Target,
                                      geneInfo = geneInfoHT,
                                      method = "geneLength")

boxplot(dataNorm, outline = FALSE)

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile",
                                  qnt.cut =  0.25)

# dataPrep_raw <- UseRaw_afterFilter(dataPrep, dataFilt)

dstyk_high
dstyk_low


v <- t(as.vector(dd2[rownames(dd2) == "ENSG00000133059",]))
plot(v)
hist(v, breaks = 50)
abline(v= 1908, col = 2)
summary(v)
v_fil <- ifelse(v > 1908, 1,0)[,1]
table(v_fil)

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,as.logical(v_fil)],
                            mat2 = dataFilt[,as.logical(v_fil)],
                            pipeline="limma",
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT", ClinicalDF = data.frame())



# DESeq -------------------------------------------------------------------

# https://rpubs.com/tiagochst/TCGAbiolinks_to_DESEq2

v <- t(as.vector(dd2[rownames(dd2) == "ENSG00000133059",]))
plot(v)
hist(v, breaks = 50)
abline(v= 1908, col = 2)
summary(v)
v_fil <- ifelse(v > 1908, 1,0)[,1]


library(DESeq2)

query.countsLUAD$dstyk <- factor(v_fil)
# data <- query.countsLUAD[,!is.na(query.countsLUAD$paper_IDH.status)]
data <- query.countsLUAD

ddsSE <- DESeqDataSet(data, design = ~ dstyk)
keep <- rowSums(counts(ddsSE)) >= 10
ddsSE <- ddsSE[keep,]
ddsSE <- DESeq(ddsSE)

resultsNames(ddsSE)
res <- results(ddsSE, name = "dstyk_1_vs_0")
dea <- as.data.frame(res)
summary(res)


# Annotation --------------------------------------------------------------
library(Homo.sapiens)
head(dea)
keytypes(Homo.sapiens)
geneid <- rownames(dea)
# geneid <- "ENSG00000066427"
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL"),
                keytype="ENSEMBL")

dim(genes)
dim(dea)

genes <- unique(genes)

tail(sort(table(genes$ENSEMBL)), 170)
table((sort(table(genes$ENSEMBL))) > 1)

# Quedarnos con un duplicado

library(dplyr)
genes %>% group_by(ENSEMBL) %>% filter(duplicated(ENSEMBL) | n()==1) -> genes2
genes2 %>% group_by(ENSEMBL) %>% filter(duplicated(ENSEMBL) | n()==1) -> genes3
genes3 %>% group_by(ENSEMBL) %>% filter(duplicated(ENSEMBL) | n()==1) -> genes4

# genes4 <- unique(genes4)
dim(genes4)
table((sort(table(genes4$ENSEMBL))) > 1)
tail(sort(table(genes4$ENSEMBL)), 15)

genes$ENSEMBL[!genes$ENSEMBL %in% rownames(dea)]
rownames(dea)[!rownames(dea) %in% genes$ENSEMBL]

dea$trans <- rownames(dea)

m1 <- merge(dea, genes4, by.x = "trans", by.y = "ENSEMBL", all.x = TRUE)

m1 %>% arrange(padj) -> m2

normalized_counts <- counts(ddsSE, normalized=TRUE)

lamp_trasn <- as.character(genes4$ENSEMBL[genes4$SYMBOL == "NQO1"])
lamp_trasn <- lamp_trasn[!is.na(lamp_trasn)]
lamp <- data.frame(normalized_counts[rownames(normalized_counts) == lamp_trasn,])
names(lamp)[1] <- "Lamp"
lamp$dstyk <- unname(v_fil)
boxplot(lamp[,1] ~ lamp$dstyk) #, ylim =c(0, 100000)

lamp %>% group_by(dstyk) %>% summarise(mean_d = mean(Lamp, na.rm = T))

save.image("20210128_all.RData")

# Survival Ibon -----------------------------------------------------------




# Survival ----------------------------------------------------------------

# Selecting only 20 genes for example
dataBRCAcomplete <- log2(normalized_counts[1:20,] + 1)

# clinical_patient_Cancer <- GDCquery_clinic("TCGA-BRCA","clinical")
clinical_patient_Cancer <- data.frame(
  bcr_patient_barcode = colnames(dataBRCAcomplete),
  vital_status = dd$vital_status,
  days_to_death = dd$days_to_death,
  days_to_last_follow_up = dd$days_to_last_follow_up
)

group1 <- TCGAquery_SampleTypes(colnames(dataBRCAcomplete), typesample = c("TP"))
group2 <- TCGAquery_SampleTypes(colnames(dataBRCAcomplete), typesample = c("NT"))

tabSurvKM <- TCGAanalyze_SurvivalKM(clinical_patient_Cancer,
                                    dataBRCAcomplete,
                                    Genelist = rownames(dataBRCAcomplete),
                                    Survresult = FALSE,
                                    p.cut = 0.4,
                                    ThreshTop = 0.67,
                                    ThreshDown = 0.33,
                                    group1 = group1, # Control group
                                    group2 = group2) # Disease group

# If the groups are not specified group1 == group2 and all samples are used
## Not run:
tabSurvKM <- TCGAanalyze_SurvivalKM(clinical_patient_Cancer,
                                    dataBRCAcomplete,
                                    Genelist = rownames(dataBRCAcomplete),
                                    Survresult = TRUE,
                                    p.cut = 0.2,
                                    ThreshTop = 0.67,
                                    ThreshDown = 0.33)

## End(Not run)


