% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/from_vcf2df.R
\name{somaticsniper_2_excel}
\alias{somaticsniper_2_excel}
\title{From vcf to excel by somatic sniper}
\usage{
somaticsniper_2_excel(
  name_vcffile,
  excel = FALSE,
  excel_file = "20200306_somaticsniper.xlsx"
)
}
\arguments{
\item{name_vcffile}{The name of the vcf file}

\item{excel}{If we want to export to excel (boolean)}

\item{excel_file}{Name of the output excel file.}
}
\value{
data frame with ordered information
}
\description{
From vcf to excel by somatic sniper
}
\examples{
\dontrun{
strelka_2_excel (name_vcffile)
}

}
