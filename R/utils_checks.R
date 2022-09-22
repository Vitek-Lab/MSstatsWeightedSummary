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
            pp_graph = createPeptideProteinGraph(feature_data)
            feature_data = addClusterMembership(feature_data, pp_graph)
        }
        feature_data
    }
}


#' Split data into a list of clusters
#' @keywords internal
getProteinsClusters = function(feature_data) {
    feature_data = feature_data[, .(Cluster, Run, ProteinName, PSM, Channel,
                                    log2IntensityNormalized)]
    split(feature_data, feature_data[["Cluster"]])
}

getAnnotation = function(feature_data) {
    unique(feature_data[, list(Run, Mixture, TechRepMixture,
                               Channel, Condition, BioReplicate)])
}
