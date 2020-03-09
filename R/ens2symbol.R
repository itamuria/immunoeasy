
#' Convert ENS to symbols
#'
#' @description Search symbols from ENG ids and calculate the length of the gene
#'
#' @param ens_ids a vector with the ensemble ids
#' @param specie Specie that we want to get information from
#' @param attributes_list attributes that we want to obtain in the last table
#' @param filter_name the source information
#'
#' @return a data framewith the symbol, chromosome, start position, end position and length
#'
#' @examples
#' ens2symbol(ens_ids = c("ENSG00000000003"))
#'
#' @export
ens2symbol <- function(ens_ids, specie = "hsapiens_gene_ensembl",
                       attributes_list = c("ensembl_gene_id","hgnc_symbol","ensembl_gene_id", "chromosome_name","start_position","end_position"),
                       filter_name = "ensembl_gene_id")
{

  # mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- biomaRt::useEnsembl("ensembl")
  mart <- biomaRt::useDataset(specie, mart)

  ens <- as.character(ens_ids)
  ens2 <- gsub('\\..+$', '', ens)

  annotLookup2 <- biomaRt::getBM(
    mart = mart,
    attributes = c("ensembl_gene_id", "hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
    filter = filter_name,
    values = ens2,
    uniqueRows = TRUE)

  annotLookup2$length <- annotLookup2$end_position - annotLookup2$start_position
  return(annotLookup2)
}
