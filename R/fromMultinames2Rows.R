#' From Multinames to rows
#'
#' @description Sometimes there are more than one names in the same field. This function helps to separate them in different rows
#'
#' @param dataset The dataset that we want to change
#' @param colu_name The name of the column that we want to take into account
#'
#' @return data frame with new rows
#' @export
#'
#' @examples
#' \dontrun {
#' varcall2HowMany (filename)
#' }
#'
from_multinames_to_rows <- function(dataset, colu_name)
{
  colu_position <- which(names(dataset) == colu_name)

  which_delete <- grep(",",dataset[,colu_name])
  dataset$delete <- 0
  dataset$delete[which_delete] <- 1

  add_1 <- NULL

  for(v in which_delete)
  {
    # print(v)
    temp <- as.character(dataset[v,colu_name])
    st <- strsplit(temp,",")
    len1 <- length(st[[1]])
    if(len1 > 1)
    {
      for(b in 1:len1)
      {
        # cufflink2$borrar[v] <- 1
        add_1 <- c(add_1, as.matrix(dataset[v,-which(names(dataset)==colu_name)]), st[[1]][b])
      }
    }
  }

  gehitu <- data.frame(matrix(add_1, ncol = ncol(dataset), byrow = T))
  # gehitu$delete <- 0
  dataset2 <- dataset[dataset$delete == 0,]

  izenak <- names(dataset)[-colu_position]
  names(gehitu) <- c(izenak, colu_name)
  gehitu$delete <- 0

  dataset3 <- rbind(dataset2, gehitu)
  return(dataset3)
}
