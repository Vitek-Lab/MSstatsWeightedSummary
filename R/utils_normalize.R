#' Normalization for TMT data
#' @param input data.table
#'
#' @return data.table
#' @export
#'
normalizeSharedPeptides = function(input) {
    Intensity = log2Intensity = NULL

    pp_match = unique(input[, list(ProteinName, PeptideSequence)])
    annotation = unique(input[, list(Run, Mixture, TechRepMixture,
                                     Channel, BioReplicate, Condition)])
    to_normalize = unique(input[, list(PeptideSequence, Charge, PSM,
                                       Run, Channel, Intensity)])
    to_normalize[, log2Intensity := log(Intensity, 2)]

    if (any(!is.na(to_normalize$Intensity) & to_normalize$Intensity < 1)) {
        to_normalize[, log2Intensity := ifelse(!is.na(Intensity) & Intensity < 1,
                                        NA, log2Intensity)]
    }
    to_normalize = .normalizePeptides(to_normalize)
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
#' @param input data.table
#' @param normalize logical, if TRUE, `input` data will be normalized
#' @return data.table
#' @keywords internal
.normalizePeptides = function(input) {
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
    # input[, log2Intensity := log2IntensityNormalized]
    input[, Intensity := 2 ^ log2IntensityNormalized]
    input[, Diff := NULL]

    input[, !(colnames(input) == "MedianLog2Int"), with = FALSE]
}
