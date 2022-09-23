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

    list(FeatureLevelData = feature_data,
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
    data.table::rbindlist(
        lapply(summary_per_cluster,
               function(cluster_summary) {
                   peptide_protein_dt = cluster_summary[["pp_dt"]]
                   final_weights_per_run = cluster_summary[["alpha_history"]]
                   data.table::rbindlist(
                       lapply(names(final_weights_per_run), function(run_id) {
                           n_iters = length(final_weights_per_run[[run_id]])
                           weights = final_weights_per_run[[run_id]][[n_iters]]
                           weights = merge(weights, peptide_protein_dt,
                                           by = c("ProteinName", "PSM"),
                                           all.x = TRUE, all.y = TRUE)
                           weights[, Weight := ifelse(is.na(Weight), 0, Weight)]
                           weights[, Run := run_id]
                           weights
                       }))
               }))
}

#' Get history of weights from all iterations
#' @inheritParams processSummarizationOutput
#' @keywords internal
getWeightsHistory = function(summary_per_cluster, save_weights_history) {
    if (save_weights_history) {
        data.table::rbindlist(
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
                                                     by = c("ProteinName", "PSM"),
                                                     all.x = T, all.y = T)
                                           })
                            iters = lapply(seq_along(iters),
                                           function(i) cbind(iters[[i]],
                                                             iter = i))
                            iters = data.table::rbindlist(iters)
                            iters[, Run := run_id]
                            iters[, Weight := ifelse(is.na(Weight), 0, Weight)]
                            iters
                        }))
                })
        )
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
                histories = data.table::rbindlist(histories_per_run)
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
