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
                   ConvergenceSummary = "data.table",
                   WeightsHistory = "dtOrNULL",
                   ConvergenceHistory = "dtOrNULL"))

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


#' Create input for MSstatsTMT::groupComparison function
#' @export
setGeneric("makeMSstatsTMTInput",
           function(weighted_summary, msstatstmt_output = NULL)
               standardGeneric("makeMSstatsTMTInput"))
#' Create input for MSstatsTMT::groupComparison function
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


#' Prepare summarizaton output
#' @inheritParams getWeightedProteinSummary
#' @param summary_per_cluster output of getClusterSummaries
#' @param annotation output of getAnnotation
#' @keywords internal
processSummarizationOutput = function(summary_per_cluster,
                                      feature_data,
                                      annotation,
                                      save_weights_history,
                                      save_convergence_history,
                                      tolerance) {
    summary = data.table::rbindlist(lapply(summary_per_cluster,
                                           function(x) x[["summary"]]))
    summary = merge(summary, annotation, by = c("Run", "Channel"))
    data.table::setnames(summary, "ProteinName", "Protein")
    weights_summary = getWeightsSummary(summary_per_cluster)
    weights_history = getWeightsHistory(summary_per_cluster,
                                        save_weights_history)

    convergence_summary = getConvergenceSummary(summary_per_cluster,
                                                tolerance)
    convergence_history = getConvergenceHistory(summary_per_cluster,
                                                tolerance,
                                                save_convergence_history)

    new("MSstatsWeightedSummary",
        FeatureLevelData = feature_data,
        ProteinLevelData = summary,
        Weights = weights_summary[order(Run, PSM)],
        ConvergenceSummary = convergence_summary,
        WeightsHistory = weights_history,
        ConvergenceHistory = convergence_history)
}

#' Get summary of final weights
#' @inheritParams processSummarizationOutput
#' @keywords internal
getWeightsSummary = function(summary_per_cluster) {
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
                                           all.x = TRUE, all.y = TRUE)
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
    weights
}

#' Get history of weights from all iterations
#' @inheritParams processSummarizationOutput
#' @keywords internal
getWeightsHistory = function(summary_per_cluster, save_weights_history) {
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
                                                     all.x = T, all.y = T)
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
        weights_history
    } else {
        NULL
    }
}

#' Get convergence summary
#' @inheritParams processSummarizationOutput
#' @keywords internal
getConvergenceSummary = function(summary_per_cluster, tolerance) {
    data.table::rbindlist(
        lapply(
            names(summary_per_cluster),
            function(cluster_summary_id) {
                cluster_summary = summary_per_cluster[[cluster_summary_id]]
                convergence_histories = cluster_summary[["convergence_history"]]
                histories_per_run = lapply(names(convergence_histories),
                                           function(run_id) {
                                               history = convergence_histories[[run_id]]
                                               list(Run = run_id,
                                                    NumIterations = length(history),
                                                    FinalDiffValue = history[length(history)],
                                                    Tolerance = tolerance,
                                                    Converged = history[length(history)] < tolerance)
                                           })
                histories = data.table::rbindlist(histories_per_run, fill = TRUE,
                                                  use.names = TRUE)
                histories[, Cluster := cluster_summary_id]
                histories[, list(Cluster, Run, NumIterations,
                                 FinalDiffValue, tolerance, Converged)]
            })
    )
}

#' Get details of convergence
#' @inheritParams processSummarizationOutput
#' @keywords internal
getConvergenceHistory = function(summary_per_cluster,
                                 tolerance,
                                 save_convergence_history) {
    if (save_convergence_history) {
        data.table::rbindlist(
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
    } else {
        NULL
    }
}
