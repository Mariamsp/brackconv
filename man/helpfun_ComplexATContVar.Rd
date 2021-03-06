% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AddComplexATContVar.R
\name{helpfun_ComplexATContVar}
\alias{helpfun_ComplexATContVar}
\title{AddComplexATContVar}
\usage{
helpfun_ComplexATContVar(df, STR_complexity_df)
}
\arguments{
\item{df}{Data frame output from \code{\link{getLUS}}}
}
\value{
Two data frames: 1) final data frame with stutter sequences and 2) data frame with noise sequences
}
\description{
Help function for adding complexity and AT-content variables to data frame
}
\section{1. Motif type as a factor variable and establish groups (do it per marker)}{

}

\section{2. STR complexity as factor variable:}{
\itemize{
\item Simple
\item Compound
\item Complex
}
}

\author{
Maria Martin Agudo <maagud.at.ous-hf.no>
}
