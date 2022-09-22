#' Normalization for TMT data with shared peptides
#'
#' @inheritParams getWeightedProteinSummary
#'
#' @return data.table
#' @export
#'
normalizeSharedPeptides = function(feature_data) {
    Intensity = log2Intensity = NULL

    pp_match = unique(feature_data[, list(ProteinName, PeptideSequence)])
    annotation = unique(feature_data[, list(Run, Mixture, TechRepMixture,
                                     Channel, BioReplicate, Condition)])
    to_normalize = unique(feature_data[, list(PeptideSequence, Charge, PSM,
                                       Run, Channel, Intensity)])
    to_normalize[, log2Intensity := log(Intensity, 2)]

    to_normalize = normalizePeptides(to_normalize)
    if (any(to_normalize$Intensity < 1 & !is.na(to_normalize$Intensity))) {
        to_normalize[, log2Intensity := ifelse(Intensity < 1 & !is.na(Intensity),
                                        NA, log2Intensity)]
        to_normalize[, Intensity := ifelse(Intensity < 1 & !is.na(Intensity),
                                    NA, Intensity)]
    }
    output = merge(to_normalize, pp_match, by = "PeptideSequence",
                   all.x = TRUE, allow.cartesian = TRUE)
    output = merge(output, annotation, by = c("Run", "Channel"))
    output
}


#' Normalization between channels (before summarization)
#' @inheritParams getWeightedProteinSummary
#' @return data.table
#' @keywords internal
normalizePeptides = function(feature_data) {
    log2Intensity = Intensity = Run = Channel = NULL
    MedianLog2Int = Diff = NULL

    input[, MedianLog2Int := median(log2Intensity, na.rm = TRUE),
          by = c("Run", "Channel")]
    median_baseline = median(
        unique(input[, list(Run, Channel, MedianLog2Int)])[, MedianLog2Int],
        na.rm = TRUE
    )
    input[, Diff := median_baseline - MedianLog2Int]
    input[, log2IntensityNormalized := log2Intensity + Diff]
    input[, Intensity := 2 ^ log2IntensityNormalized]
    input[, Diff := NULL]

    input[, !(colnames(input) == "MedianLog2Int"), with = FALSE]
}
