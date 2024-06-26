% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TCGA_dataDownload.R
\name{getTCGA_CNV.data}
\alias{getTCGA_CNV.data}
\title{getTCGA_CNV.data}
\usage{
getTCGA_CNV.data(
  project,
  save = FALSE,
  folder = ".",
  data.type = "Gene Level Copy Number"
)
}
\arguments{
\item{project}{A list of valid project (see list with TCGAbiolinks:::getGDCprojects()$project_id)]}

\item{save}{TRUE or FALSE}

\item{folder}{The path specifies the location of the folder where the data is saved when the 'save' parameter is set to TRUE.}

\item{data.type}{"Gene Level Copy Number" or "Gene Level Copy Number Scores"}
}
\description{
getTCGA_CNV.data
}
