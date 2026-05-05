#' Check if data is label-free or TMT
#' @keywords internal
checkExperimentType = function(feature_data) {
    if (is.element("Channel", colnames(feature_data))) {
        "TMT"
    } else {
        "LF"
    }
}

#' Check if data is in MSstatsTMT format
#' @keywords internal
checkDataCorrectness = function(feature_data, experiment_type) {
    `:=` = Cluster = Run = IsotopeLabelType = FragmentIon = ProductCharge = NULL
    Fraction = PSM = PeptideSequence = PrecursorCharge = NULL

    if (experiment_type == "TMT") {
        required_columns = c("ProteinName", "PeptideSequence", "Charge",
                             "PSM", "Channel", "Intensity", "Run",
                             "Condition", "BioReplicate",
                             "Mixture", "TechRepMixture")
        optional_columns = c("Cluster", "log2IntensityNormalized")
        all_columns = c(required_columns, optional_columns)

        if (!all(required_columns %in% colnames(feature_data))) {
            stop(paste("Missing columns in input.",
                       "Please verify that data is in the MSstatsTMT format"))
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
    } else {
        required_columns = c("ProteinName", "PeptideSequence", "PrecursorCharge",
                             "Intensity", "Run",
                             "Condition", "BioReplicate")
        optional_columns = c("Cluster", "log2IntensityNormalized",
                             "FragmentIon", "ProductCharge", "Fraction",
                             "IsotopeLabelType", "IsUnique", "PSM")
        all_columns = c(required_columns, optional_columns)

        if (!all(required_columns %in% colnames(feature_data))) {
            stop(paste("Missing columns in input.",
                       "Please verify that data is in the MSstats format"))
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
            if (!is.element("IsotopeLabelType", colnames(feature_data))) {
                feature_data[, IsotopeLabelType := "Light"]
            }
            if (!is.element("FragmentIon", colnames(feature_data))) {
                feature_data[, FragmentIon := NA_character_]
            }
            if (!is.element("ProductCharge", colnames(feature_data))) {
                feature_data[, ProductCharge := NA_character_]
            }
            if (!is.element("Fraction", colnames(feature_data))) {
                feature_data[, Fraction := "1"]
            }
            if (!is.element("PSM", colnames(feature_data))) {
                feature_data[, PSM := paste(PeptideSequence, PrecursorCharge,
                                            FragmentIon, ProductCharge,
                                            sep = "_")]
            }

            # feature_data = imputeTMT(feature_data)
            feature_data

        }
    }
}

#' Split data into a list of clusters
#' @inheritParams getWeightedProteinSummary
#' @keywords internal
getProteinsClusters = function(feature_data) {
    Cluster = Run = ProteinName = PSM = Channel = log2IntensityNormalized = NULL
    feature_data = feature_data[, list(Cluster, Run, ProteinName, PSM, Channel,
                                    log2IntensityNormalized)]
    split(feature_data, feature_data[["Cluster"]])
}

#' Make annotation
#' @keywords internal
getAnnotation = function(feature_data, experiment_type) {
    Run = Mixture = TechRepMixture = Channel = Condition = BioReplicate = NULL
    Fraction = IsotopeLabelType = NULL
    if (experiment_type == "TMT") {
        unique(feature_data[, list(Run, Mixture, TechRepMixture,
                                   Channel, Condition, BioReplicate)])
    } else {
        unique(feature_data[, list(Run, Fraction, IsotopeLabelType,
                                   Condition, BioReplicate)])
    }
}

#' @keywords internal
getLFDataPortion = function(feature_data, experiment_type) {
    PSM = PeptideSequence = PrecursorCharge = FragmentIon = ProductCharge = NULL
    Run = Fraction = IsotopeLabelType = NULL

    if (experiment_type == "TMT") {
        NULL
    } else {
        unique(feature_data[, list(PSM, PeptideSequence, PrecursorCharge,
                                   FragmentIon, ProductCharge, Run, Fraction,
                                   IsotopeLabelType)])
    }
}

#' @keywords internal
reshapeLFData = function(feature_data, experiment_type) {
    Cluster = ProteinName = PeptideSequence = PrecursorCharge = PSM = NULL
    Condition = BioReplicate = Run = Intensity = log2IntensityNormalized = NULL

    if (experiment_type == "TMT") {
        feature_data
    } else {
        feature_data[, list(
            Cluster,
            ProteinName,
            PeptideSequence,
            Charge = PrecursorCharge,
            PSM,
            Mixture = "1",
            TechRepMixture = "1",
            Run = "1",
            Condition,
            BioReplicate,
            Channel = Run,
            Intensity,
            log2IntensityNormalized)]
    }
}
