#' @importClassesFrom data.table data.table
setClassUnion("dtOrNULL", c("data.table", "NULL"))

#' Output of weighted summarization
#'
#' @slot FeatureLevelData feature-level (input) data
#' @slot ProteinLevelData protein-level (summarized) output data
#' @slot Weights a table of final peptide-protein Weights
#' @slot ConvergenceSummary table with information about convergence for each Cluster and Run
#' @slot WeightsHistory optional data.table of Weights from all iterations of fitting algorithm
#' @slot ConvergenceHistory optional data.table with sums of absolute values of differences between Weights from consecutive iteration
#'
setClass("MSstatsWeightedSummary",
         slots = c(FeatureLevelData = "data.table",
                   ProteinLevelData = "data.table",
                   Weights = "data.table",
                   FittedProfiles = "data.table",
                   ConvergenceSummary = "data.table",
                   FinalCriterionValues = "data.table",
                   WeightsHistory = "dtOrNULL",
                   ConvergenceHistory = "dtOrNULL",
                   ExperimentType = "character"))

#' Extract feature-level data from MSstatsWeightedSummary object
#' @export
setGeneric("featureData",
           function(weighted_summary, proteins = NULL) standardGeneric("featureData"))
#' Extract feature-level data from MSstatsWeightedSummary object
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @param proteins optional character vector of proteins to extract. If NULL, all proteins
#' will be returned
#' @return data.table
setMethod("featureData", "MSstatsWeightedSummary",
          function(weighted_summary, proteins = NULL) {
              feature_level =  weighted_summary@FeatureLevelData
              if (!is.null(proteins)) {
                  feature_level = feature_level[ProteinName %in% proteins]
              }
              feature_level
          })


#' Extract protein-level data from MSstatsWeightedSummary object
#' @export
setGeneric("proteinData",
           function(weighted_summary, proteins = NULL) standardGeneric("proteinData"))
#' Extract protein-level data from MSstatsWeightedSummary object
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @param proteins optional character vector of proteins to extract. If NULL, all proteins
#' will be returned
#' @return data.table
setMethod("proteinData", "MSstatsWeightedSummary",
          function(weighted_summary, proteins = NULL) {
              protein_level =  weighted_summary@ProteinLevelData
              if (!is.null(proteins)) {
                  protein_level = protein_level[ProteinName %in% proteins]
              }
              protein_level
          })

#' Extract weights data from MSstatsWeightedSummary object
#' @export
setGeneric("featureWeights",
           function(weighted_summary, proteins = NULL, shared_only = TRUE)
               standardGeneric("featureWeights"))
#' Extract weights data from MSstatsWeightedSummary object
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @param proteins optional character vector of proteins to extract. If NULL, all proteins
#' will be returned
#' @param shared_only logical, output data.table will only include shared peptides
#' @return data.table
setMethod("featureWeights", "MSstatsWeightedSummary",
          function(weighted_summary, proteins = NULL, shared_only = TRUE) {
              weights = weighted_summary@Weights
              if (!is.null(proteins)) {
                  weights = weights[ProteinName %in% proteins]
              }
              if (shared_only) {
                  weights = weights[!(IsUnique)]
              }
              weights
          })

#' Extract convergence information from MSstatsWeightedSummary object
#' @export
setGeneric("convergenceSummary",
           function(weighted_summary) standardGeneric("convergenceSummary"))
#' Extract convergence information from MSstatsWeightedSummary object
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @return data.table
setMethod("convergenceSummary", "MSstatsWeightedSummary",
          function(weighted_summary) {
              weighted_summary@ConvergenceSummary
          })

#' Extract weights history from MSstatsWeightedSummary object
#' @export
setGeneric("weightsHistory",
           function(weighted_summary, shared_only = TRUE) standardGeneric("weightsHistory"))
#' Extract weights history from MSstatsWeightedSummary object
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @param shared_only logical, output data.table will only include shared peptides
#' @return data.table
setMethod("weightsHistory", "MSstatsWeightedSummary",
          function(weighted_summary, shared_only = TRUE) {
              weights_history = weighted_summary@WeightsHistory
              if (shared_only) {
                  weights_history = weights_history[!(IsUnique)]
              }
              weights_history
          })

#' Extract convergence history from MSstatsWeightedSummary object
#' @export
setGeneric("convergenceHistory",
           function(weighted_summary) standardGeneric("convergenceHistory"))
#' Extract convergence history from MSstatsWeightedSummary object
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @return data.table
setMethod("convergenceHistory", "MSstatsWeightedSummary",
          function(weighted_summary) {
              weighted_summary@ConvergenceHistory
          })

#' Extract cluster information from MSstatsWeightedSummary object
#' @export
setGeneric("proteinClusters",
           function(weighted_summary) standardGeneric("proteinClusters"))
#' Extract cluster information from MSstatsWeightedSummary object
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @return data.table
setMethod("proteinClusters", "MSstatsWeightedSummary",
          function(weighted_summary) {
              feature_data = weighted_summary@FeatureLevelData
              cluster_data = unique(feature_data[, .(Run, Cluster, ProteinName)])
              cluster_data
          })

#' Extract fitted PSM-level profiles from MSstatsWeightedSummary object
#' @export
setGeneric("fittedProfiles",
           function(weighted_summary) standardGeneric("fittedProfiles"))
#' Extract cluster information from MSstatsWeightedSummary object
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @return data.table
setMethod("fittedProfiles", "MSstatsWeightedSummary",
          function(weighted_summary) {
              fitted_profiles = weighted_summary@FittedProfiles
              fitted_profiles
          })


#' Plot summary for a given cluster or proteins set
#' @export
#' @return ggplot2 object
setGeneric("plotSummary",
           function(weighted_summary, cluster = NULL,
                    proteins = NULL, channel_order = NULL)
               standardGeneric("plotSummary"))
#' Plot summary for a given cluster or proteins set
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @param cluster optional ID of a cluster to plot. Either `cluster` or `proteins` must be provided
#' @param proteins optional vector of proteins to plot. Either `cluster` or `proteins` must be provided
#' @param channel_order optional vector of ordered channel IDs. If provided,
#' x-axis of the plot will follow this order
#' @return ggplot2 object
#' @import ggplot2
setMethod("plotSummary", "MSstatsWeightedSummary",
          function(weighted_summary, cluster = NULL,
                   proteins = NULL, channel_order = NULL) {
              if (is.null(cluster) & is.null(proteins)) {
                  stop("cluster or proteins must be provided")
              } else {
                  if (is.null(cluster)) {
                      feature_data = weighted_summary@FeatureLevelData[ProteinName %in% proteins]
                  } else {
                      feature_data = weighted_summary@FeatureLevelData[Cluster == cluster]
                      proteins = feature_data[, unique(ProteinName)]
                  }
                  protein_data = weighted_summary@ProteinLevelData[Protein %in% proteins]
                  data.table::setnames(protein_data, "Protein", "ProteinName")

                  if (!is.null(channel_order)) {
                      feature_data[, Channel := factor(Channel, channel_order,
                                                       ordered = TRUE)]
                      protein_data[, Channel := factor(Channel, channel_order,
                                                       ordered = TRUE)]
                  }

                  ggplot(protein_data, aes(x = Channel, y = Abundance,
                                           group = ProteinName, color = ProteinName)) +
                      geom_line(aes(x = Channel, y = log2IntensityNormalized,
                                    group = PSM), data = feature_data,
                                inherit.aes = FALSE) +
                      geom_line(size = 1.5) +
                      facet_grid(Run ~ ProteinName) +
                      theme_bw() +
                      theme(legend.position = "bottom")
              }
          })


#' Create input for MSstatsTMT::groupComparisonTMT function
#' @export
setGeneric("makeMSstatsTMTInput",
           function(weighted_summary, msstatstmt_output = NULL)
               standardGeneric("makeMSstatsTMTInput"))
#' Create input for MSstatsTMT::groupComparisonTMT function
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @param msstatstmt_output optional output of MSstatsTMT::proteinSummarization function
#' @return list
setMethod("makeMSstatsTMTInput", "MSstatsWeightedSummary",
          function(weighted_summary, msstatstmt_output = NULL) {
              feature_data = weighted_summary@FeatureLevelData
              protein_data = weighted_summary@ProteinLevelData

              if (is.null(msstatstmt_output)) {
                  list(FeatureLevelData = feature_data,
                       ProteinLevelData = protein_data)
              } else {
                  mstmt_feature = msstatstmt_output[["FeatureLevelData"]]
                  mstmt_protein = msstatstmt_output[["ProteinLevelData"]]
                  list(FeatureLevelData = rbind(feature_data,
                                                mstmt_feature,
                                                use.names = TRUE, fill = TRUE),
                       ProteinLevelData = rbind(protein_data,
                                                mstmt_protein,
                                                use.names = TRUE, fill = TRUE))
              }
          })

#' Create input for MSstats::groupComparison function
#' @export
setGeneric("makeMSstatsInput",
           function(weighted_summary, msstats_output = NULL)
               standardGeneric("makeMSstatsInput"))
#' Create input for MSstats::groupComparison function
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @param msstatstmt_output optional output of MSstats::dataProcess function
#' @return list
setMethod("makeMSstatsInput", "MSstatsWeightedSummary",
          function(weighted_summary, msstats_output = NULL) {
              feature_data = data.table::copy(weighted_summary@FeatureLevelData)
              protein_data = data.table::copy(weighted_summary@ProteinLevelData)

              if (!is.null(msstats_output)) {
                  mst_feature = msstats_output[["FeatureLevelData"]]
                  mst_protein = msstats_output[["ProteinLevelData"]]
                  feature_data_conv = rbind(feature_data,
                                       mst_feature,
                                       use.names = TRUE, fill = TRUE)
                  protein_data = rbind(protein_data,
                                       mst_protein,
                                       use.names = TRUE, fill = TRUE)
              }
              if (weighted_summary@ExperimentType == "LF") {
                  feature_data_conv = feature_data[, list(
                      PROTEIN = ProteinName,
                      PEPTIDE = paste(PeptideSequence, Charge, sep = "_"),
                      TRANSITION = paste(FragmentIon, ProductCharge, sep = "_"),
                      FEATURE = PSM,
                      LABEL = IsotopeLabelType,
                      GROUP = Condition,
                      RUN = Run,
                      SUBJECT = BioReplicate,
                      FRACTION = Fraction,
                      originalRUN = Run,
                      censored = FALSE,
                      INTENSITY = Intensity,
                      ABUNDANCE = log2IntensityNormalized,
                      newABUNDANCE = log2IntensityNormalized,
                                                     predicted = NA_real_,
                                                     remove = FALSE)]
                  cols = c("RUN", "Protein", "LogIntensities", "originalRUN",
                           "GROUP", "SUBJECT", "more50missing", "NumMeasuredFeature")
                  setnames(protein_data,
                           c("Run", "Abundance", "Condition", "BioReplicate"),
                           c("originalRUN", "LogIntensities", "GROUP", "SUBJECT"))
                  protein_data[, RUN := originalRUN]
                  protein_data[, more50missing := FALSE]
                  num_features = feature_data[!is.na(log2IntensityNormalized),
                                              list(NumMeasuredFeature = data.table::uniqueN(PSM)),
                                              by = c("ProteinName", "Run")]
                  group_meas = feature_data[!is.na(log2IntensityNormalized),
                                            list(TotalGroupMeasurements = data.table::uniqueN(PSM)),
                                            by = c("ProteinName", "Condition")]
                  protein_data = merge(protein_data, num_features,
                                       by.x = c("Protein", "originalRUN"),
                                       by.y = c("ProteinName", "Run"),
                                       all.x = TRUE, all.y = TRUE, sort = FALSE)
                  protein_data = merge(protein_data, group_meas,
                                       by.x = c("Protein", "GROUP"),
                                       by.y = c("ProteinName", "Condition"),
                                       sort = FALSE)
                  protein_data[, MissingPercentage := 0.0]
                  protein_data[, NumImputedFeature := 0L]
              } else {
                  feature_data_conv = feature_data
              }
              list(FeatureLevelData = feature_data_conv,
                   ProteinLevelData = protein_data)
          })


#' Prepare summarizaton output
#' @inheritParams getWeightedProteinSummary
#' @param summary_per_cluster output of getClusterSummaries
#' @param annotation output of getAnnotation
#' @keywords internal
processSummarizationOutput = function(summary_per_cluster,
                                      feature_data,
                                      lf_data,
                                      annotation,
                                      save_weights_history,
                                      save_convergence_history,
                                      tolerance,
                                      experiment_type) {
    summary = combineSummaries(summary_per_cluster, annotation,
                               experiment_type)

    fitted_profiles = getAllFittedProfiles(summary_per_cluster,
                                           experiment_type)

    weights_summary = getWeightsSummary(summary_per_cluster,
                                        experiment_type)
    weights_history = getWeightsHistory(summary_per_cluster,
                                        save_weights_history)
    criteria = getFinalCriteria(summary_per_cluster, experiment_type)

    convergence_summary = getConvergenceSummary(summary_per_cluster,
                                                tolerance,
                                                experiment_type)
    convergence_history = getConvergenceHistory(summary_per_cluster,
                                                tolerance,
                                                save_convergence_history,
                                                experiment_type)

    final_feature_data = getFinalFeatureData(feature_data,
                                             lf_data,
                                             experiment_type)

    new("MSstatsWeightedSummary",
        FeatureLevelData = final_feature_data,
        ProteinLevelData = summary,
        Weights = weights_summary[order(Run, PSM)],
        FittedProfiles = fitted_profiles,
        ConvergenceSummary = convergence_summary,
        FinalCriterionValues = criteria,
        WeightsHistory = weights_history,
        ConvergenceHistory = convergence_history,
        ExperimentType = experiment_type)
}

#' @keywords internal
combineSummaries = function(summary_per_cluster, annotation,
                            experiment_type) {
    summary = data.table::rbindlist(lapply(summary_per_cluster,
                                           function(x) x[["summary"]]))
    if (experiment_type == "TMT") {
        summary = merge(summary, annotation, by = c("Run", "Channel"), sort = FALSE)
    } else {
        summary[, Run := NULL]
        data.table::setnames(summary, "Channel", "Run")
        summary = merge(summary, annotation, by = "Run", sort = FALSE)
    }
    data.table::setnames(summary, "ProteinName", "Protein")
    summary
}

#' Get summary of final weights
#' @inheritParams processSummarizationOutput
#' @keywords internal
getWeightsSummary = function(summary_per_cluster, experiment_type) {
    weights = data.table::rbindlist(
        lapply(summary_per_cluster,
               function(cluster_summary) {
                   peptide_protein_dt = cluster_summary[["pp_dt"]]
                   final_weights_per_run = cluster_summary[["alpha_history"]]
                   data.table::rbindlist(
                       lapply(names(final_weights_per_run), function(run_id) {
                           n_iters = length(final_weights_per_run[[run_id]])
                           weights = final_weights_per_run[[run_id]][[n_iters]]
                           weights[, Run := run_id]
                           weights = merge(weights, peptide_protein_dt,
                                           by = c("ProteinName", "PSM", "Run"),
                                           all.x = TRUE, all.y = TRUE, sort = FALSE)
                           weights[, Weight := ifelse(is.na(Weight), 0, Weight)]
                           weights[, Total := sum(Weight),
                                   by = c("ProteinName", "Run")]
                           weights = weights[Total > 0]
                           weights[, Total := NULL]
                           weights
                       }), fill = TRUE, use.names = TRUE)
               }), fill = TRUE, use.names = TRUE)
    weights[, IsUnique := data.table::uniqueN(ProteinName) == 1,
            by = c("PSM", "Run")]
    if (experiment_type == "LF") {
        weights[, Run := NA_character_]
    }
    weights
}

#' Get history of weights from all iterations
#' @inheritParams processSummarizationOutput
#' @keywords internal
getWeightsHistory = function(summary_per_cluster, save_weights_history,
                             experiment_type) {
    if (save_weights_history) {
        weights_history = data.table::rbindlist(
            lapply(
                summary_per_cluster,
                function(cluster_summary) {
                    peptide_protein_dt = cluster_summary[["pp_dt"]]
                    weight_history = cluster_summary[["alpha_history"]]
                    data.table::rbindlist(
                        lapply(names(weight_history), function(run_id) {
                            iters = weight_history[[run_id]]
                            iters = lapply(iters,
                                           function(x) {
                                               merge(x,
                                                     peptide_protein_dt,
                                                     by = c("ProteinName", "PSM", "Run"),
                                                     all.x = T, all.y = T, sort = FALSE)
                                           })
                            iters = lapply(seq_along(iters),
                                           function(i) cbind(iters[[i]],
                                                             iter = i))
                            iters = data.table::rbindlist(iters, fill = TRUE, use.names = TRUE)
                            iters[, Weight := ifelse(is.na(Weight), 0, Weight)]
                            iters
                        }), fill = TRUE, use.names = TRUE)
                }), fill = TRUE, use.names = TRUE
        )
        weights_history[, IsUnique := data.table::uniqueN(ProteinName) == 1,
                        by = c("PSM", "Run")]
        if (experiment_type == "LF") {
            weights_history[, Run := NA_character_]
        }
        weights_history
    } else {
        NULL
    }
}

#' Get convergence summary
#' @inheritParams processSummarizationOutput
#' @keywords internal
getConvergenceSummary = function(summary_per_cluster, tolerance,
                                 experiment_type) {
    conv_summary =data.table::rbindlist(
        lapply(
            names(summary_per_cluster),
            function(cluster_summary_id) {
                histories_per_run = lapply(names(summary_per_cluster[[cluster_summary_id]][["convergence_history"]]),
                                           function(run_id) {
                                               history = summary_per_cluster[[cluster_summary_id]][["convergence_history"]][[run_id]]
                                               list(Run = run_id, NumIterations = length(summary_per_cluster[[cluster_summary_id]][["alpha_history"]][[run_id]]),
                                                    FinalDiffValue = ifelse(length(history[length(history)]) == 0, 0.0, history[length(history)]),
                                                    Tolerance = tolerance,
                                                    Converged = ifelse(length(history[length(history)]) == 0, TRUE, history[length(history)] <
                                                                           tolerance))
                                           })
                histories = data.table::rbindlist(histories_per_run, fill = TRUE,
                                                  use.names = TRUE)
                histories[, Cluster := cluster_summary_id]
                histories[, list(Cluster, Run, NumIterations,
                                 FinalDiffValue, tolerance, Converged)]
            })
    )
    if (experiment_type == "LF") {
        conv_summary[, Run := NA_character_]
    }
    conv_summary
}

#' Get details of convergence
#' @inheritParams processSummarizationOutput
#' @keywords internal
getConvergenceHistory = function(summary_per_cluster,
                                 tolerance,
                                 save_convergence_history,
                                 experiment_type) {
    if (save_convergence_history) {
        conv_history = data.table::rbindlist(
            lapply(
                names(summary_per_cluster),
                function(cluster_summary_id) {
                    cluster_summary = summary_per_cluster[[cluster_summary_id]]
                    convergence_histories = cluster_summary[["convergence_history"]]
                    histories_per_run = lapply(names(convergence_histories),
                                               function(run_id) {
                                                   history = convergence_histories[[run_id]]
                                                   list(Run = run_id,
                                                        Iter = seq_along(history),
                                                        DiffValue = history)
                                               })
                    histories = data.table::rbindlist(histories_per_run)
                    histories[, Cluster := cluster_summary_id]
                    histories[, NumIterations := max(Iter),
                              by = "Run"]
                    histories[, Converged := min(DiffValue) < tolerance,
                              by = "Run"]
                    histories[, list(Cluster, Run, Iter, DiffValue,
                                     NumIterations, Converged)]
                })
        )
        if (experiment_type == "LF") {
            conv_history[, Run := NA_character_]
        }
        conv_history
    } else {
        NULL
    }
}

#' @keywords internal
getFinalCriteria = function(summary_per_cluster, experiment_type) {
    criteria = data.table::rbindlist(lapply(names(summary_per_cluster), function(cluster_id) {
        run_summaries = summary_per_cluster[[cluster_id]]
        cbind(Cluster = cluster_id,
              run_summaries[["final_criterion_values"]])
    }))
    data.table::setnames(criteria, "Criterion", "FinalCriterion")
    if (experiment_type == "LF") {
        criteria[, Run := NA_character_]
    }
    criteria
}

#' @keywords internal
getAllFittedProfiles = function(summary_per_cluster, experiment_type) {
    data.table::rbindlist(lapply(names(summary_per_cluster), function(cluster_id) {
        run_summaries = summary_per_cluster[[cluster_id]]
        fitted_profiles = run_summaries[["estimated_profiles"]]
        pp_dt = run_summaries[["pp_dt"]]
        fitted_profiles = merge(fitted_profiles, pp_dt,
                                by = c("Run", "PSM"),
                                all.x = TRUE, all.y = TRUE, sort = FALSE,
                                allow.cartesian = TRUE)
        fitted_profiles[, Cluster := cluster_id]
        fitted_profiles = fitted_profiles[, list(Cluster, ProteinName, PSM, Run,
                                                 Channel, log2IntensityNormalized, Predicted)]
        if (experiment_type == "LF") {
            fitted_profiles[, Run := NULL]
            data.table::setnames(fitted_profiles, "Channel", "Run")
        }
        fitted_profiles
    }))
}

#' @keywords internal
getFinalFeatureData = function(feature_data, lf_data, experiment_type) {
    if (experiment_type == "LF") {
        feature_data = merge(feature_data, lf_data,
                             by.x = c("Channel", "PeptideSequence", "Charge", "PSM"),
                             by.y = c("Run", "PeptideSequence", "PrecursorCharge", "PSM"))
        feature_data[, Run := NULL]
        data.table::setnames(feature_data, "Channel", "Run")
    }
    feature_data
}
