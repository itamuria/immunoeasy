
# https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
# https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

setwd("G:/Mi unidad/NASIR/NASIR_counts")
source("functions.r")


pat <- openxlsx::read.xlsx("../MAIN_FOLDER_NASIR/MAIN_TABLE_v002.xlsx",2)
load("202101_counts_subread_together.RData")

names(temp)
head(pat)

h <- temp[,c(1,15:44)]

library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)

file_s <- dir("Raw_Counts/")
file2 <- grep("subread_tumor_rna_counts", file_s, value = T)
file3 <- grep("summary", file2, value = T, invert = T)
file4 <- grep("RNAP1", file3, value = T, invert = T)
file5 <- paste0("Raw_Counts/",file4)

k <- readDGE(file5, columns=c(1,7), header = T, skip = 1)
class(k)
dim(k)


# Clustering --------------------------------------------------------------


scal_count <- scale(k$counts)

seeds_df_sc <- as.data.frame(scale(t(k$counts)))
# summary(seeds_df_sc)
dist_mat <- dist(seeds_df_sc, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
par(mfrow=c(1,1))
plot(hclust_avg)
cut_avg <- cutree(hclust_avg, k = 4)
plot(hclust_avg)
rect.hclust(hclust_avg , k = 4, border = 2:6)
abline(h = 3, col = 'red')
suppressPackageStartupMessages(library(dendextend))
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, h = 3)
plot(avg_col_dend)

suppressPackageStartupMessages(library(dplyr))
seeds_df_cl <- mutate(seeds_df_sc, cluster = cut_avg)
count(seeds_df_cl,cluster)

library(pheatmap)
# pheatmap(scale(k$counts),main = "pheatmap default", scale = "row",show_rownames = F)


# Creating object ---------------------------------------------------------



samplenames <- substring(colnames(k), 19, nchar(colnames(k)))
samplenames <- gsub("_subread_tumor_rna_counts","",samplenames)
samplenames

colnames(k) <- samplenames
group <- as.factor(c("Pro1",
                     "Pro0","Pro0","Pro1","Pro0","Pro1","Pro1","Pro1","Pro1",
                     "Pro0", "Pro0","Pro0",
                     "Pro1","Pro1","Pro0","Pro0","Pro0","Pro1","Pro1","Pro0","Pro0","Pro1","Pro1","Pro1","Pro1",
                     "Pro1","Pro1","Pro1", "Pro1","Pro0"))
k$samples$group <- group
# lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
# k$samples$lane <- lane
k$samples

# head(read.delim(file5[1], header = TRUE, skip = 1))
# read.delim(file5[1], skip = 1)

trans <- read.table("transcripts_symbol.txt")
names(trans) <- c("Chr","Pos1","Pos2","ENSG","VV","ENST","V","trans_type","VVV","Symbol","sasf","gene_type","kkk")
trans <- trans[,c(1:4,6,8,10,12)]


duplicated(c(1,1,1,3))
k$genes <- trans

cpm <- cpm(k)
lcpm <- cpm(k, log=TRUE)

L <- mean(k$samples$lib.size) * 1e-6
M <- median(k$samples$lib.size) * 1e-6
c(L, M)

table(rowSums(k$counts==0)==9)

keep.exprs <- filterByExpr(k, group=group)
k <- k[keep.exprs,, keep.lib.sizes=FALSE]
dim(k)


# calcNormFactors(k)

lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(k)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(k, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

k <- calcNormFactors(k, method = "TMM")
k$samples$norm.factors

x2 <- k
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5


par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)
x2$samples$norm.factors

lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")


lcpm <- cpm(k, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Progression")
plotMDS(lcpm, labels=rownames(k$samples), dim=c(1,2))
title(main="B. Samples")


# DEG ---------------------------------------------------------------------


design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

contr.matrix <- makeContrasts(
  Pro_01 = Pro0 - Pro1,
  levels = colnames(design))
contr.matrix

par(mfrow=c(1,2))
v <- voom(k, design, plot=TRUE)
v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

summary(decideTests(efit))

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

top <- topTreat(tfit, coef=1, n=Inf)
openxlsx::write.xlsx(top, file = "top_progresion.xlsx")

# Cluster Normalizado -----------------------------------------------------

scal_count <- scale(t(cpm(k, log=TRUE)))

seeds_df_sc <- as.data.frame(scal_count)
# summary(seeds_df_sc)
dist_mat <- dist(seeds_df_sc, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
par(mfrow=c(1,1))
plot(hclust_avg)
cut_avg <- cutree(hclust_avg, k = 4)
plot(hclust_avg)
rect.hclust(hclust_avg , k = 4, border = 2:6)
abline(h = 3, col = 'red')
suppressPackageStartupMessages(library(dendextend))
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, h = 3)
plot(avg_col_dend)

pheatmap(cpm(k, log=TRUE),main = "pheatmap default", scale = "row",show_rownames = F)

