#' Check if data is in MSstatsTMT format
#' @keywords internal
checkDataCorrectness = function(feature_data) {
    required_columns = c("ProteinName", "PeptideSequence", "Charge",
                         "PSM", "Channel", "Intensity", "Run",
                         "Condition", "BioReplicate",
                         "Mixture", "TechRepMixture")
    optional_columns = c("Cluster", "log2IntensityNormalized")
    all_columns = c(required_columns, optional_columns)

    if (!all(required_columns %in% colnames(feature_data))) {
        stop(paste("Missing columns in input.",
                   "Please verify that data is in MSstatsTMT format"))
    } else {
        feature_data = data.table::as.data.table(feature_data)
        if (!is.element("log2IntensityNormalized", colnames(feature_data))) {
            feature_data = normalizeSharedPeptides(feature_data)
        }

        all_columns = intersect(all_columns, colnames(feature_data))
        feature_data = feature_data[, all_columns, with = FALSE]

        if (!is.element("Cluster", colnames(feature_data))) {
            by_run = split(feature_data, feature_data[["Run"]])
            by_run = lapply(by_run, function(x) {
                pp_graph = createPeptideProteinGraph(x)
                x = addClusterMembership(x, pp_graph)
                x
            })
            feature_data = data.table::rbindlist(by_run)
            feature_data[, Cluster := paste(Cluster, Run, sep = "__")]
        }
        # feature_data = imputeTMT(feature_data)
        feature_data
    }
}


#' Split data into a list of clusters
#' @inheritParams getWeightedProteinSummary
#' @keywords internal
getProteinsClusters = function(feature_data) {
    feature_data = feature_data[, .(Cluster, Run, ProteinName, PSM, Channel,
                                    log2IntensityNormalized)]
    split(feature_data, feature_data[["Cluster"]])
}

#' Make annotation
#' @keywords internal
getAnnotation = function(feature_data) {
    unique(feature_data[, list(Run, Mixture, TechRepMixture,
                               Channel, Condition, BioReplicate)])
}
