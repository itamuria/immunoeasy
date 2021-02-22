library("DESeq2")
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(officer)
library(EnhancedVolcano)

val_contr <- NULL
contraste = "Progresion"
val_padj = 0.1

sample <- data.frame(main[,c("ID","TUMOR_RNA", contraste)])

sample[,contraste] <- factor(sample[,contraste])
sele1 <- which(!is.na(sample[,contraste]))
length(sele1)
sample$sele1 <- 888
sample$sele1[sele1] <- "SI"
sele2 <- which(!is.na(sample$TUMOR_RNA))
length(sele2)
sample$sele2 <- 888
sample$sele2[sele2] <- "SI"
sele <- ifelse(sample$sele1 == "SI" & sample$sele2 == "SI", TRUE, FALSE)

tum_rna <- sample$TUMOR_RNA[!is.na(sample$TUMOR_RNA)]
se <- se2[,which(colnames(se2) %in% tum_rna)]
contr <- sample$TUMOR_RNA[sele]
se3 <- se[,which(colnames(se) %in% contr)]


# Solo tumores ------------------------------------------------------------
var_formula <- as.formula(paste0("~ ", contraste))
dds <- DESeqDataSetFromMatrix(countData = se3, colData = sample[sele,], design = var_formula)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)


se3_froga <- se3[apply(se3, 1, sum) > 5,]
boxplot(log(se3+1))
boxplot(log(se3_froga+1))
summary(normalized_counts)
boxplot(log(normalized_counts+1))

# write.table(normalized_counts, file=paste0("SandraIbon/",contraste,"_normalized_counts.txt"),
#             sep="\t", quote=F, col.names=NA)

# Diff expression
dds <- DESeq(dds)
res <- results(dds)
# head(results(dds, tidy=TRUE)) #let's look at the results table
summary(res)

res <- res[order(res$padj),]
head(res)

res_frame <- data.frame(trans = rownames(res), data.frame(res))

symbol <- mapIds(org.Hs.eg.db,
                 keys=gsub("\\..*","",rownames(res_frame)),
                 column="SYMBOL",
                 keytype="ENSEMBL",
                 multiVals="first")
symbol <- unname(ifelse(is.na(symbol), names(symbol), symbol))
res_frame2 <- data.frame(sym = data.frame(symbol), res_frame)


# openxlsx::write.xlsx(res_frame2, file = paste0("ResultsCounts2/",contraste,"_DEG_with_Full_.xlsx"))


# Volcano -----------------------------------
mm <- -log10(min(res_frame2$padj, na.rm = TRUE))

volc1 <- EnhancedVolcano(res_frame2,
                         lab = res_frame2$symbol,
                         x = 'log2FoldChange',
                         y = 'padj',
                         # xlim = c(-8, 8),
                         ylim = c(0,mm),
                         title = paste0(contraste, '_with_padj01'),
                         pCutoff = 10e-2,
                         FCcutoff = 1,
                         pointSize = 3.0,
                         labSize = 3.0)


res_val3 <- res_frame2[res_frame2$padj < val_padj,]
res_val3 <- res_val3[!is.na(res_val3$symbol),]
openxlsx::write.xlsx(res_val3, file = paste0("ResultsCounts2/",contraste,"_DEG_with_Filt_.xlsx"))

# ClusterProfiler

res_val3 <- res_val3[order(res_val3$log2FoldChange, decreasing = TRUE),]
d <- res_val3[,c("symbol","log2FoldChange")]
geneList <- d[,2]
names(geneList) <- as.character(d[,1])
