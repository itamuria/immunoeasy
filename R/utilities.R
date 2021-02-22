#' Extract type of extension
#'
#' @param filename Filename that need to be extracted
#'
#' @return last three letters (extension)
#' @export
#'
#' @examples
#' \dontrun{
#' filter_format_word ("testing.doc")
#' }
#'
filter_format_word <- function(filename)
{
  lenword <- nchar(filename)
  initialpoint <- lenword - 2
  return(substr(filename,initialpoint,lenword))
}


#' Select those elements of a vector with the specified extension
#'
#' @param word_vector Vector with several names
#' @param formatw Extension to be filtered
#'
#' @return Select the elements that fill the criteria
#' @export
#'
#' @examples
#' \dontrun{
#' select_files_format (c("testing.vcf","bat.vcf","bibi.cod"), formatw = "vcf")
#' }
#'
select_files_format <- function(word_vector, formatw = "vcf")
{
  new_vector <- NULL
  for(h in 1:length(word_vector))
  {
    if(filter_format_word(word_vector[h]) == formatw) new_vector <- c(new_vector, word_vector[h])
  }
  return(new_vector)
}


#' Extract day, hour or day and hour of the moment
#'
#' @param time1 day, hour or day_hour
#'
#' @return A character with the information of the moment in different formats to include in names of files
#' @export
#'
#' @examples
#'
#' \dontrun{
#' day_hour ("day")
#' day_hour ("hour")
#' day_hour ("day_hour")
#' }
day_hour <- function(time1) {
  if (time1 == "day") {
    day <- paste(strsplit(as.character(substring(Sys.time(), 1, 10)), "-")[[1]], collapse = "")
    fetxia <-  day

  } else if (time1 == "hour") {
    hour <- paste0(substring(Sys.time(), 12, 13), "h", substring(Sys.time(), 15, 16), "m", substring(Sys.time(), 18,19), "s")
    fetxia <- hour

  } else if (time1 == "day_hour") {
    day <- paste(strsplit(as.character(substring(Sys.time(), 1, 10)), "-")[[1]], collapse = "")
    hour <- paste0(substring(Sys.time(), 12, 13), "h", substring(Sys.time(), 15, 16), "m", substring(Sys.time(), 18, 19), "s")
    fetxia <- paste0(day, "_", hour)

  } else {
    print("Choose the right time, please")
  }
  return(fetxia)
}

#' Save pheatmap to pdf
#'
#' @param x pheatmap object
#' @param filename filename
#' @param width width
#' @param height height
#'
#' @return Save the pheatmap in a pdf
#' @export
#'
#' @examples
#'
#' \dontrun{
#' save_pheatmap_pdf (x, filename, width=7, height=7)
#' }
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
