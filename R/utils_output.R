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
