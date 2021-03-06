#' Merge all the counters together
#'
#' @param subread_name Data frame with normalized subread cuartiles
#' @param cuff_name Data frame with normalized cufflink cuartiles
#' @param htseq_name Data frame with normalized htseq cuartiles
#' @param quant_name Data frame with normalized quant cuartiles
#' @param ngenes Selection of genes that we want to filter
#' @param mfl_number Mean length
#' @param export_excel_name Name of the final excel name
#' @param save_intermediate If we want to save intermediate files
#' @param dif_cuartiles If we want to calculate the differences between cuartiles
#'
#' @return data frame with all the information merged
#' @export
#'
#' @examples
#' \dontrun{
#' four_counter2summary (subread_name,cuff_name,quant_name,htseq_name,
#'                                ngenes = NULL, mfl_number = mfl_number,
#'                                export_excel_name = "four_together.xlsx",
#'                                save_intermediate = FALSE)
#' }
#'
four_counter2summary <- function(semi_subread,semi_cuff,semi_quant,semi_htseq,
                                ngenes = NULL,
                                export_excel_name = "four_together.xlsx",
                                save_final = FALSE, dif_cuartiles = FALSE)
{
  names(semi_cuff) <- c("Gene_id","Symbol","cuff_locus","cuff_FPKM","cuff_Cuartiles")
  semi_cuffb10 <- from_multinames_to_rows (dataset = semi_cuff, colu_name = "Symbol")
  semi_cuffb10 <- semi_cuffb10[,-which(names(semi_cuffb10)=="delete")]

  semi_htseq2 <- semi_htseq[,c(1,6,5,2,4)]
  semi_htseq2$ht_locus <- paste0(semi_htseq$chromosome_name,":",semi_htseq$start_position,"-",semi_htseq$end_position)
  names(semi_htseq2) <- c("Gene_id","Symbol","ht_Counts","ht_FPKM","ht_Cuartiles","ht_locus")
  semi_htseq2 <- semi_htseq2[which(semi_htseq2$ht_locus%in%grep("^C|^G",semi_htseq2$ht_locus, value=TRUE, invert = TRUE)),]

  semi_quant2 <- semi_quant[,c(1,5,2,4)]
  names(semi_quant2) <- c("Symbol","Quan_Counts","Quan_FPKM","Quan_Cuartiles")

  semi_subread2 <- semi_subread[,c(1,6,5,2,4)]
  semi_subread2$ht_locus <- paste0(semi_subread$chromosome_name,":",semi_subread$start_position,"-",semi_subread$end_position)
  names(semi_subread2) <- c("Gene_id","Symbol","SR_Counts","SR_FPKM","SR_Cuartiles","SR_locus")
  semi_subread2 <- semi_subread2[which(semi_subread2$SR_locus%in%grep("^C|^G",semi_subread2$SR_locus, value=TRUE, invert = TRUE)),]

  # merging

  m1_logico <- merge(semi_subread2[semi_subread2$Symbol != "",], semi_quant2, by = "Symbol")
  m2_logico <- merge(semi_cuffb10, semi_htseq2[semi_htseq2$Symbol != "",], by = "Symbol", all.y = TRUE)
  m2_logico$Symbol <- as.character(m2_logico$Symbol)
  m3_logico <- merge(m1_logico, m2_logico, by = "Symbol", all = TRUE)

  # filtrar con nuestros datos o no.

  mfinal <- m3_logico[m3_logico$Symbol%in%ngenes,]

  mfinal2 <- mfinal[,c("Symbol","Gene_id","Gene_id.x","SR_locus","ht_locus","cuff_locus","SR_Counts","Quan_Counts","ht_Counts",
                       "SR_FPKM","Quan_FPKM","cuff_FPKM","ht_FPKM","SR_Cuartiles","Quan_Cuartiles","cuff_Cuartiles","ht_Cuartiles")]

  mfinal2$SR_Cuartiles <- as.numeric(as.character(mfinal2$SR_Cuartiles))
  mfinal2$Quan_Cuartiles <- as.numeric(as.character(mfinal2$Quan_Cuartiles))
  mfinal2$cuff_Cuartiles <- as.numeric(as.character(mfinal2$cuff_Cuartiles))
  mfinal2$ht_Cuartiles <- as.numeric(as.character(mfinal2$ht_Cuartiles))

  if(dif_cuartiles)
  {
    mfinal2$dif_SR_Quan <- as.numeric(as.character(mfinal2$SR_Cuartiles)) - as.numeric(as.character(mfinal2$Quan_Cuartiles))
    mfinal2$dif_SR_ht <- as.numeric(as.character(mfinal2$SR_Cuartiles)) - as.numeric(as.character(mfinal2$ht_Cuartiles))
    mfinal2$dif_SR_cuff <- as.numeric(as.character(mfinal2$SR_Cuartiles)) - as.numeric(as.character(mfinal2$cuff_Cuartiles))
    mfinal2$dif_Quan_ht <- as.numeric(as.character(mfinal2$Quan_Cuartiles)) - as.numeric(as.character(mfinal2$ht_Cuartiles))
    mfinal2$dif_Quan_cuff <- as.numeric(as.character(mfinal2$Quan_Cuartiles)) - as.numeric(as.character(mfinal2$cuff_Cuartiles))
    mfinal2$dif_ht_cuff <- as.numeric(as.character(mfinal2$ht_Cuartiles)) - as.numeric(as.character(mfinal2$cuff_Cuartiles))
  }

  if(save_final)
  {
    openxlsx::write.xlsx(mfinal2, file=export_excel_name)
  }

  return(mfinal2)

}
