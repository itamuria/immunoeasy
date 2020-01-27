#' Fpkm to quartiles
#'
#' @param fpkm_values
#'
#' @return a vector with quartiles
#' @export
#'
#' @examples
#' fpkm2cuartiles_cuff(rnorm(100))
#'
fpkm2cuartiles <- function(fpkm_values)
{
  qt <- quantile(fpkm_values[fpkm_values!=0], probs=0:4/4)
  q1 <- as.numeric(qt[1])
  q2 <- as.numeric(qt[2])
  q3 <- as.numeric(qt[3])
  q4 <- as.numeric(qt[4])

  qtal <- ifelse(fpkm_values <= q2, 1,
                 ifelse(fpkm_values <= q3 & fpkm_values > q2, 2,
                        ifelse(fpkm_values <= q4 & fpkm_values > q3, 3,
                               ifelse(fpkm_values > q4, 4, NA))))
  return(qtal)
}
