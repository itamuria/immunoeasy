#' Merge all the counters together
#'
#' @param subread_raw_file Data frame with normalized subread cuartiles
#' @param cufflink_raw_file Data frame with normalized cufflink cuartiles
#' @param htseq_raw_file Data frame with normalized htseq cuartiles
#' @param quant_raw_file Data frame with normalized quant cuartiles
#'
#' @return data frame with all the information merged
#' @export
#'
#' @examples
#' \dontrun {
#' four_counter2summary (subread_name,cuff_name,quant_name,htseq_name,
#'                                ngenes, mfl_number = mfl_number, export_excel_name = "four_together.xlsx", save_intermediate = FALSE)
#' }
#'
four_counter2summary <- function(subread_name,cuff_name,quant_name,htseq_name,
                                ngenes, mfl_number = mfl_number, export_excel_name = "four_together.xlsx", save_intermediate = FALSE)
{
    semi_subread <- counts2fpkm_subread (subread_name, mfl_num = c(mfl_number))
    semi_cuff <- counts2fpkm_cuff (subread_name)
    semi_cuant <- counts2fpkm_quant (quant_name, mfl_num = c(mfl_number))
    semi_htseq <- counts2fpkm_htseq (htseq_name, mfl_num = c(mfl_number))

  head(semi_cuff)
  semi_cuffb <- semi_cuff[,c(1,5,7,10,14)]
  names(semi_cuffb) <- c("Gene_id","Symbol","cuff_locus","cuff_FPKM","cuff_Cuartiles")
  head(semi_cuffb)
  semi_cuffb10 <- from_multinames_to_rows (dataset = semi_cuffb, colu_name = "Symbol")
  semi_cuffb10 <- semi_cuffb10[,-which(names(semi_cuffb10)=="delete")]


  head(semi_cuff2)
  semi_cuff2b <- semi_cuff2[,c(1,5,7,10,14)]
  names(semi_cuff2b) <- c("Gene_id","Symbol","cuff2_locus","cuff2_FPKM","cuff2_Cuartiles")
  head(semi_cuff2b)
  semi_cuff2b10 <- from_multinames_to_rows (dataset = semi_cuff2b, colu_name = "Symbol")
  semi_cuff2b10 <- semi_cuff2b10[,-which(names(semi_cuff2b10)=="delete")]


  head(semi_htseq)
  semi_htseq2 <- semi_htseq[,c(1,7,2,4,6)]
  semi_htseq2$ht_locus <- paste0(semi_htseq$chromosome_name,":",semi_htseq$start_position,"-",semi_htseq$end_position)
  names(semi_htseq2) <- c("Gene_id","Symbol","ht_FPKM","ht_Cuartiles","ht_Counts","ht_locus")
  head(semi_htseq2)
  semi_htseq2 <- semi_htseq2[which(semi_htseq2$ht_locus%in%grep("^C|^G",semi_htseq2$ht_locus, value=TRUE, invert = TRUE)),]

  head(semi_quant)
  names(semi_quant) <- c("Symbol","Quan_FPKM","Quan_Cuartiles","Quan_Counts")
  head(semi_quant)

  head(semi_subread)
  semi_subread2 <- semi_subread[,c(1,7,2,4,6)]
  semi_subread2$ht_locus <- paste0(semi_subread$chromosome_name,":",semi_subread$start_position,"-",semi_subread$end_position)
  names(semi_subread2) <- c("Gene_id","Symbol","SR_FPKM","SR_Cuartiles","SR_Counts","SR_locus")
  head(semi_subread2)
  semi_subread2 <- semi_subread2[which(semi_subread2$SR_locus%in%grep("^C|^G",semi_subread2$SR_locus, value=TRUE, invert = TRUE)),]

  # merging

  m1_total <- merge(semi_subread2, semi_quant, by = "Symbol", all = TRUE)
  # quitando los vacioes en symbol
  m1_logico <- merge(semi_subread2[semi_subread2$Symbol!="",], semi_quant, by = "Symbol")

  semi_cuff_c <- data.frame(semi_cuff2b10,semi_cuffb10$cuff_Cuartiles[semi_cuffb10$Symbol!="-"])
  head(semi_cuff_c)
  names(semi_cuff_c)[6] <- "cuff_Cuartiles"
  semi_cuff_c$cuff2_locus <- as.character(semi_cuff_c$cuff2_locus)
  semi_cuff_c <- semi_cuff_c[which(semi_cuff_c$cuff2_locus%in%grep("^H|^G|^C",semi_cuff_c$cuff2_locus, value=TRUE, invert = TRUE)),]

  # m2_total <- merge(semi_cuff_c, semi_htseq2, by = "Symbol", all = TRUE)
  # quitando los vacioes en symbol
  # table(semi_cuffb$Symbol!="-")
  m2_logico <- merge(semi_cuff_c, semi_htseq2, by = "Symbol", all.y = TRUE)

  # m3_total <- merge(m1_total,m2_total, by = "Symbol", all = TRUE)
  m3_logico <- merge(m1_logico,m2_logico, by = "Symbol", all = TRUE)

  # openxlsx::write.xlsx(m2_total, file="subread_quant_merged.xlsx")
  # openxlsx::write.xlsx(m1_total, file="cuff_htseq_merged.xlsx")

  # openxlsx::write.xlsx(m3_logico, file="four_together.xlsx")

  # filtrar con nuestros datos



  mfinal <- m3_logico[m3_logico$Symbol%in%ngenes,]
  ngenes[!ngenes%in%m3_logico$Symbol]

  head(mfinal)
  mfinal2 <- mfinal[,c(1,2,15,6,11,19,5,9,18,3,7,12,16,4,8,13,14,17)]
  head(mfinal2)


  mfinal2$SR_Cuartiles <- as.numeric(as.character(mfinal2$SR_Cuartiles))
  mfinal2$Quan_Cuartiles <- as.numeric(as.character(mfinal2$Quan_Cuartiles))
  mfinal2$cuff_Cuartiles <- as.numeric(as.character(mfinal2$cuff_Cuartiles))
  mfinal2$cuff2_Cuartiles <- as.numeric(as.character(mfinal2$cuff2_Cuartiles))
  mfinal2$ht_Cuartiles <- as.numeric(as.character(mfinal2$ht_Cuartiles))

  mfinal2$dif_SR_Quan <- as.numeric(as.character(mfinal2$SR_Cuartiles)) - as.numeric(as.character(mfinal2$Quan_Cuartiles))
  mfinal2$dif_SR_ht <- as.numeric(as.character(mfinal2$SR_Cuartiles)) - as.numeric(as.character(mfinal2$ht_Cuartiles))
  mfinal2$dif_SR_cuff <- as.numeric(as.character(mfinal2$SR_Cuartiles)) - as.numeric(as.character(mfinal2$cuff_Cuartiles))
  mfinal2$dif_SR_cuff2 <- as.numeric(as.character(mfinal2$SR_Cuartiles)) - as.numeric(as.character(mfinal2$cuff2_Cuartiles))

  mfinal2$dif_Quan_ht <- as.numeric(as.character(mfinal2$Quan_Cuartiles)) - as.numeric(as.character(mfinal2$ht_Cuartiles))
  mfinal2$dif_Quan_cuff <- as.numeric(as.character(mfinal2$Quan_Cuartiles)) - as.numeric(as.character(mfinal2$cuff_Cuartiles))
  mfinal2$dif_Quan_cuff2 <- as.numeric(as.character(mfinal2$Quan_Cuartiles)) - as.numeric(as.character(mfinal2$cuff2_Cuartiles))
  mfinal2$dif_ht_cuff <- as.numeric(as.character(mfinal2$ht_Cuartiles)) - as.numeric(as.character(mfinal2$cuff_Cuartiles))
  mfinal2$dif_ht_cuff2 <- as.numeric(as.character(mfinal2$ht_Cuartiles)) - as.numeric(as.character(mfinal2$cuff2_Cuartiles))
  mfinal2$dif_cuff_cuff2 <- as.numeric(as.character(mfinal2$cuff_Cuartiles)) - as.numeric(as.character(mfinal2$cuff2_Cuartiles))

  openxlsx::write.xlsx(mfinal2, file=export_excel_name)
  return(list(mfinal2=mfinal2, semi_cuff_c = semi_cuff_c, semi_htseq2 = semi_htseq2,
              semi_subread2 = semi_subread2, semi_quant = semi_quant))

}
