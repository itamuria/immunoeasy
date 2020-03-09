#' From fpkm to quartiles with Cufflink files
#'
#' @param filename The name of the file with Cufflink counts
#' @param previous_clean If we should remove the gene_names == "-" before FPKM calculus
#'
#' @return data frame with fpkm and quartiles
#' @export
#'
#' @examples
#' \dontrun{
#' counts2fpkm_cuff (filename)
#' }
#'
counts2fpkm_cuff <- function(filename, previous_clean = FALSE)
{
  # Load data  ------------------------------------------------------------
  dat <- utils::read.table(filename, header = T, sep = "\t")
  dat <- dat[,c("gene_id", "gene_short_name", "locus", "FPKM")]

  if(previous_clean){
    dat <- dat[dat$gene_short_name != "-",]
  }
  dat$Cuartiles <- fpkm2cuartiles(dat$FPKM)

  return(dat)
}
