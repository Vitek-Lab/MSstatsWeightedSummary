#' @description
#' MSstatsWeightedSummary implements a statistical model for joint analysis of proteins that share peptides based on mass spectrometry data.
#' The package extends existing MSstats and MSstatsTMT workflows to include shared peptides in input data for protein summarization.
#' The summarization workflow is implemented in the \code{\link{getWeightedProteinSummary}} function. For example, check the package vignette
#' `browseVignettes(package = "MSstatsWeightedSummary.")`.
#' For the statistical details and background, please refer to the publication \link{https://doi.org/10.1093/bioinformatics/btaf046}.
#' @keywords internal
"_PACKAGE"

#' Simulated data set illustrating quantification with shared peptides
#'
#' Quantification data in the MSstatsTMT format with characteristics similar to the
#' BRD case study presented in the original publication.
#'
#' @format A data frame (data.table) with 168 rows and 14 variables:
#' \describe{
#'   \item{ProteinName}{label of protein matching feature's AA sequence}
#'   \item{PeptideSequence}{AA sequence of the feature}
#'   \item{Charge}{charge of the feature)}
#'   \item{PSM}{PSM label - concatenation of PeptideSequence and Charge columns}
#'   \item{Run}{label of MS run (constant)}
#'   \item{Mixture}{label of a TMT mixture (constant)}
#'   \item{TechRepMixture}{label of technical replicate mixture (constant)}
#'   \item{Channel}{label of a TMT channel}
#'   \item{BioReplicate}{label of biological replicate (group comparison design)}
#'   \item{Condition}{label of a group}
#'   \item{Log2Intensity}{log-intensity values}
#'   \item{Intensity}{intensity values}
#'   \item{IsUnique}{logical column indicating if a given peptide is unique to its parent protein}
#'   \item{log2IntensityNormalized}{normalized log-intensities of features}
#' }
#'
"simulated_dataset"