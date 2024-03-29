% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/counts2fpkm_subread.R
\name{counts2fpkm_subread}
\alias{counts2fpkm_subread}
\title{From counts to fpkm and quartiles with subread}
\usage{
counts2fpkm_subread(filename, mfl_num = c(mfl_number))
}
\arguments{
\item{filename}{The name of the file with subread counts}

\item{mfl_num}{meanFragmentLength. A numeric vector with mean fragment lengths, which can be calculated using ’CollectInsertSizeMetrics(Picard)’ tool. The length of items should be as the same of columns in read count matrix.}
}
\value{
data frame with fpkm and quartiles
}
\description{
From counts to fpkm and quartiles with subread
}
\examples{
\dontrun{
counts2fpkm_subread (filename, mfl_num = c(mfl_number))
}

}
