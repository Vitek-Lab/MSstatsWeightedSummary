% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_weights.R
\name{getPeptideProteinWeights}
\alias{getPeptideProteinWeights}
\title{Get PSM-protein weights for summarization with shared peptides}
\usage{
getPeptideProteinWeights(
  feature_data,
  norm = "p_norm",
  norm_parameter = 1,
  weights_mode = "contributions",
  weights_penalty = FALSE,
  weights_penalty_param = 0.1
)
}
\arguments{
\item{feature_data}{data.table in MSstatsTMT format. See also the Details section}

\item{norm}{"p_norm" or "Huber"}

\item{norm_parameter}{p for norm=="p_norm", M for norm=="Huber"}

\item{weights_mode}{"contributions" for "sum to one" and "non-negative" conditions,
"probabilities" for only "non-negative" condition.}

\item{weights_penalty}{if TRUE, weights will be penalized for deviations from equal value
for all proteins matching to a given PSM}

\item{weights_penalty_param}{penalty parameter}
}
\value{
data.table
}
\description{
Get PSM-protein weights for summarization with shared peptides
}
