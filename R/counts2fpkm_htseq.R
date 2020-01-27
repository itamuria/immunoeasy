#' From counts to fpkm and quartiles with htseq files
#'
#' @param filename The name of the file with htseq counts
#' @param mfl_num meanFragmentLength. A numeric vector with mean fragment lengths, which can be calculated using ’CollectInsertSizeMetrics(Picard)’ tool. The length of items should be as the same of columns in read count matrix.
#'
#' @return data frame with fpkm and quartiles
#' @export
#'
#' @examples
#' \dontrun {
#' counts2fpkm_htseq (filename, mfl_num = c(mfl_number))
#' }
#'
counts2fpkm_htseq <- function(filename, mfl_num = c(mfl_number))
{

  # Load data  ------------------------------------------------------------
  dat <- read.table(filename)
  names(dat) <- c("GeneID","Counts")

  # Obtain symbols and length  ------------------------------------------------------------

  library("biomaRt")
  annotLookup2 <- ens2symbol(dat$GeneID)

  # ENS ID generalizing   ------------------------------------------------------------

  dat$GeneID <- gsub('\\..+$', '', dat$GeneID)

  # Remove duplicate agregating as sum of copies
  masde1 <- names(table(dat$GeneID)[table(dat$GeneID)>1])
  temp3 <- dat[dat$GeneID%in%masde1,]
  temp3b <- aggregate(temp3$Counts, by = list(temp3$GeneID), sum)
  unique_name <- unique(temp3$GeneID)
  dat2 <- dat[-which(dat$GeneID%in%unique_name),]
  names(temp3b) <- names(dat)
  dat <- rbind(dat2, temp3b)

  rownames(dat) <- dat[,1]

  gene_names2 <- merge(dat, annotLookup2, by.x = "GeneID", by.y = "ensembl_gene_id")


  # Prepare for FPKM ------------------------------------------------------------

  library(countToFPKM)

  # repetir si solo hay una columna de conteos
  mfl <- c(rep(mfl_num,2))

  # annot3$ensembl_gene_id[annot3$ensembl_gene_id%in%rownames(total4b)]
  nam_delete <- names(table(gene_names2$ensembl_gene_id)[table(gene_names2$ensembl_gene_id)>1])
  gene_names2 <- unique(gene_names2)

  gene_names3 <- gene_names2[!gene_names2$ensembl_gene_id%in%nam_delete,]

  temp6 <- gene_names2[gene_names2$ensembl_gene_id%in%nam_delete,]
  temp6$GeneID <- paste0(temp6$GeneID,c("a","b"))
  gene_names4 <- rbind(gene_names3,temp6)

  df_calc <- data.frame(gene_names4$Counts,gene_names4$Counts)
  rownames(df_calc) <- gene_names4$GeneID
  colnames(df_calc) <- c("Sample_1","Sample_2")

  ## fpkm ------------------------------------------------------------

  fpkm_matrix <- countToFPKM::fpkm (df_calc, gene_names4$length, mfl)
  fpkm_matrix2 <- data.frame(fpkm_matrix)
  fpkm_matrix2$Gene_id <- rownames(fpkm_matrix)
  fpkm_matrix2$Cuartiles <- fpkm2cuartiles(fpkm_matrix2$Sample_1)

  fpkm_matrix3 <- merge(fpkm_matrix2, gene_names4, by.x = "Gene_id", by.y = "GeneID")
}
