#' From counts to fpkm and quartiles with Quant3 files
#'
#' @param filename The name of the file with Quant3 counts
#' @param mfl_num meanFragmentLength. A numeric vector with mean fragment lengths, which can be calculated using ’CollectInsertSizeMetrics(Picard)’ tool. The length of items should be as the same of columns in read count matrix.
#'
#' @return data frame with fpkm and quartiles
#' @export
#'
#' @examples
#' \dontrun{
#' counts2fpkm_quant (filename, mfl_num = c(mfl_number))
#' }
#'
counts2fpkm_quant <- function(filename, mfl_num = c(mfl_number))
{
  # Load data  ------------------------------------------------------------

  dat <- read.table(filename, header = TRUE)
  dat <- data.frame(rownames(dat),dat[,1])
  names(dat) <- c("GeneID","Counts")

  # Update with biomart  ------------------------------------------------------------

  annotLookup2 <- ens2symbol(dat$GeneID,attributes_list=c("ensembl_gene_id","hgnc_symbol","ensembl_gene_id", "chromosome_name","start_position","end_position"),
    filter_name="hgnc_symbol")

  # ENS ID generalizing   ------------------------------------------------------------

  rownames(dat) <- dat[,1]
  gene_names2 <- merge(dat, annotLookup2, by.x = "GeneID", by.y = "hgnc_symbol")
  gene_names2 <- gene_names2[gene_names2$chromosome_name%in%c(1:22,"X","Y"),]

  # Prepare for FPKM ------------------------------------------------------------

  library(countToFPKM)

  # repetir si solo hay una columna de conteos
  mfl <- c(rep(mfl_num,2))
  delete_str <- c("__alignment_not_unique","__ambiguous","__no_feature","__not_aligned","__too_low_aQual")
  gene_names2$GeneID <- as.character(gene_names2$GeneID)
  gene_names2 <- gene_names2[!gene_names2$GeneID%in%delete_str,]

  nam_delete <- names(table(gene_names2$GeneID)[table(gene_names2$GeneID)>1])

  gene_names2 <- unique(gene_names2)
  if( length(nam_delete) > 0 )
  {
    gene_names3 <- gene_names2[!gene_names2$GeneID%in%nam_delete,]

    temp6 <- gene_names2[gene_names2$GeneID%in%nam_delete,]
    temp6$GeneID <- paste0(temp6$GeneID,c("a","b"))
    gene_names4 <- rbind(gene_names3,temp6)
  } else {
    gene_names4 <- gene_names3
  }

  df_calc <- data.frame(gene_names4$Counts,gene_names4$Counts)
  rownames(df_calc) <- gene_names4$GeneID
  colnames(df_calc) <- c("Sample_1","Sample_2")

  ## fpkm ------------------------------------------------------------

  fpkm_matrix <- countToFPKM::fpkm (df_calc, gene_names4$length, mfl)
  fpkm_matrix2 <- data.frame(fpkm_matrix)
  fpkm_matrix2$Gene_id <- rownames(fpkm_matrix)
  fpkm_matrix2$Cuartiles <- fpkm2cuartiles(fpkm_matrix2$Sample_1)

  fpkm_matrix3 <- merge(fpkm_matrix2, gene_names4, by.x = "Gene_id", by.y = "GeneID")
  return(fpkm_matrix3)
}
