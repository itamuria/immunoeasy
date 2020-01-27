#' From fpkm to quartiles with Cufflink files
#'
#' @param filename The name of the file with Cufflink counts
#'
#' @return data frame with fpkm and quartiles
#' @export
#'
#' @examples
#' \dontrun {
#' counts2fpkm_cuff (filename)
#' }
#'
counts2fpkm_cuff <- function(filename)
{

  # Load data  ------------------------------------------------------------
  dat <- read.table(filename, header = T, sep = "\t")
  dat <- dat[,c("gene_id", "gene_short_name", "locus", "FPKM")]
  return(dat)
}
