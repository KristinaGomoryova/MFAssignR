% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RecalFunctions.R
\name{FilterAndMerge}
\alias{FilterAndMerge}
\title{FilterAndMerge function}
\usage{
FilterAndMerge(Recal, column_indices, column_name)
}
\arguments{
\item{Recal}{The input data frame.}

\item{column_indices}{A vector of column indices to be selected.}

\item{column_name}{The name of the column used for filtering.}
}
\value{
A modified data frame with selected columns and filtered rows.
}
\description{
This function takes a data frame 'Recal', a vector of column indices 'column_indices',
and a column name 'column_name'. It performs the following steps:
\enumerate{
\item Selects specified columns from 'Recal'.
\item Filters rows with complete cases for the specified column.
\item Merges the original 'Recal' data frame with the filtered columns.
}
}
