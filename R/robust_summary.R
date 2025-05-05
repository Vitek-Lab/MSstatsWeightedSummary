#' Get robust protein-level summary based on unique and shared peptides
#'
#' @param feature_data data.table in MSstatsTMT format. See also the Details section
#' @param norm "p_norm" or "Huber"
#' @param norm_parameter p for norm=="p_norm", M for norm=="Huber"
#' @param weights_mode "contributions" for "sum to one" and "non-negative" conditions,
#' "probabilities" for only "non-negative" condition.
#' @param tolerance tolerance to indicate weights convergence
#' @param max_iter maximum number of iteration of the procedure
#' @param initial_summary "unique", "flat" or "flat unique"
#' @param weights_penalty if TRUE, weights will be penalized for deviations from equal value
#' for all proteins matching to a given PSM
#' @param weights_penalty_param penalty parameter
#' @param save_weights_history logical, if TRUE, weights from all iterations will
#' be returned
#' @param save_convergence_history logical, if TRUE, all differences between consecutive
#' weights estimator from all iterations will be returned
#'
#' @return list of data frames with summary and other information. See the Details section for more information
#'
#' @details
#' 1. Input format: this function takes as input data in MSstatsTMT format,
#' which is a data frame with columns ProteinName, PeptideSequence, Charge,
#' PSM (equal to PeptideSequence and Charge separated by an underscore),
#' Channel, Intensity, Run and annotation columns: BioReplicate, Condition,
#' Mixture, and TechRepMixture. Additionally, we use two columns:
#' log2IntensityNormalized and Cluster. The first column stores log-transformed
#' normalized intensities (which can be obtained with \code{\link{normalizeSharedPeptides}} function).
#' If this column is not provided, data will be normalized before summarization.
#' The second column stores information about connected sub-graphs of the
#' peptide-protein graph. This column can be added with \code{\link{addClusterMembership}} function or
#' omitted. In the second case, this information will be added before summarization.
#'
#' 2. Output format: an S4 object of class "MSstatsWeightedSummary" which consists of the following items:
#' \itemize{
#' \item{FeatureLevelData:}{feature-level (input) data}
#' \item{ProteinLevelData:}{protein-level (summarized) output data}
#' \item{Weights:}{a table of final peptide-protein Weights}
#' \item{ConvergenceSummary:}{table with information about convergence for each Cluster and Run}
#' \item{WeightsHistory:}{optional data.table of Weights from all iterations of fitting algorithm}
#' \item{ConvergenceHistory:}{optional data.table with sums of absolute values of differences between Weights from consecutive iteration}
#' }
#' Elements of this object can be accessed with functions
#' \code{\link{featureData}}, \code{\link{proteinData}},
#' \code{\link{featurWeights}}, \code{\link{convergenceSummary}},
#' \code{\link{weightsHistory}}, \code{\link{convergenceHistory}}
#'
#' For statistical details about the method, please consult the vignette.
#'
#' @export
#'
getWeightedProteinSummary = function(feature_data,
                                     norm = "p_norm", norm_parameter = 1,
                                     weights_mode = "contributions",
                                     tolerance = 1e-1, max_iter = 10,
                                     initial_summary = "unique",
                                     weights_penalty = FALSE,
                                     weights_penalty_param = 0.1,
                                     save_weights_history = FALSE,
                                     save_convergence_history = FALSE) {
    feature_data = checkDataCorrectness(feature_data)
    annotation = getAnnotation(feature_data)
    cluster_input = getProteinsClusters(feature_data)

    summary_per_cluster = getClusterSummaries(cluster_input,
                                              norm, norm_parameter,
                                              weights_mode,
                                              tolerance, max_iter,
                                              initial_summary,
                                              weights_penalty,
                                              weights_penalty_param)
    summaries = processSummarizationOutput(summary_per_cluster,
                                           feature_data,
                                           annotation,
                                           save_weights_history,
                                           save_convergence_history,
                                           tolerance)
    summaries
}

#' Calculate weighted summaries for a list of clusters
#' @param cluster_input list of protein clusters
#' @inheritParams getWeightedProteinSummary
#' @keywords internal
getClusterSummaries = function(cluster_input,
                               norm, norm_parameter,
                               weights_mode,
                               tolerance, max_iter,
                               initial_summary,
                               weights_penalty,
                               weights_penalty_param) {
    lapply(
        cluster_input,
        function(single_cluster) {
            input_by_run = split(single_cluster, single_cluster[, Run])

            output_by_run = lapply(input_by_run, function(x) {
                peptide_protein_dt = unique(x[, .(ProteinName, PSM, Run)])
                getWeightedSummarySingleRun(x, peptide_protein_dt,
                                            norm, norm_parameter,
                                            weights_mode,
                                            tolerance, max_iter,
                                            initial_summary,
                                            weights_penalty,
                                            weights_penalty_param)
            })

            summarized_output = data.table::rbindlist(
                lapply(output_by_run,
                       function(x) x[["summary"]]))
            alphas_list = lapply(output_by_run,
                                 function(x) x[["alpha_history"]])
            alpha_diffs = lapply(
                output_by_run, function(x) x[["convergence_history"]])

            list(summary = summarized_output,
                 pp_dt = peptide_protein_dt,
                 alpha_history = alphas_list,
                 convergence_history = alpha_diffs)
        })
}

#' Robust summary for a single run
#' @param peptide_protein_dt table of peptide-protein matches
#' @inheritParams getWeightedProteinSummary
#' @keywords internal
getWeightedSummarySingleRun = function(feature_data, peptide_protein_dt,
                                       norm, norm_parameter, weights_mode,
                                       tolerance, max_iter, initial_summary,
                                       weights_penalty, weights_penalty_param) {
    weights_list = vector("list", max_iter)
    input_loop = feature_data[, .(Run, ProteinName, PSM,
                                  Channel,
                                  log2IntensityNormalized)]
    num_weights = data.table::uniqueN(peptide_protein_dt)
    current_weights = rep(1, num_weights)
    previous_weights = rep(0, num_weights)

    initial_summary = getInitialSummary(input_loop,
                                        norm, norm_parameter, initial_summary)
    input_loop = merge(input_loop,
                       initial_summary[, list(Run, ProteinName, Channel,
                                              Abundance, CenteredAbundance)],
                       by = c('Run', "ProteinName", "Channel"))
    iter = 1
    weights_diffs = vector("numeric", max_iter)
    while (sum(abs(current_weights - previous_weights)) > tolerance) {
        previous_weights = current_weights
        weights = getPeptideProteinWeights(input_loop, norm,
                                           norm_parameter, weights_mode,
                                           weights_penalty, weights_penalty_param)
        weights_list[[iter]] = weights

        input_loop = merge(input_loop,
                           unique(weights[, .(PSM, ProteinName)]),
                           by = c("ProteinName", "PSM"))
        # in case some weights were eliminated

        new_abundances = summarizeProteinsClusterSingleRun(input_loop,
                                                           weights,
                                                           norm, norm_parameter,
                                                           TRUE)
        input_loop = merge(input_loop[, list(ProteinName, PSM, Channel,
                                             log2IntensityNormalized)],
                           new_abundances, by = c("ProteinName", "Channel"))
        current_weights = getCurrentWeights(weights, peptide_protein_dt)

        if (iter >= max_iter) {
            current_weights = previous_weights
        } else {
            weights_diffs[iter] = sum(abs(current_weights - previous_weights))
            iter = iter + 1
        }
    }

    summarized_output = new_abundances[, list(ProteinName, Run, Channel, Abundance)]
    weights_diffs = weights_diffs[weights_diffs > 0]
    weights_list = weights_list[sapply(weights_list, function(x) !is.null(x))]
    list(summary = summarized_output,
         alpha_history = weights_list,
         convergence_history = weights_diffs)
}

getInitialSummary = function(input_loop,
                             norm, norm_parameter,
                             initial_summary) {
    if (initial_summary == "unique") {
        initial_weights = unique(input_loop[, .(ProteinName, PSM, Weight = 1)])
        summarizeProteinsClusterSingleRun(input_loop,
                                          initial_weights,
                                          norm, norm_parameter,
                                          use_shared = FALSE)
    } else if (initial_summary == "flat") {
        means = input_loop[, .(Abundance = mean(log2IntensityNormalized,
                                                na.rm = TRUE)),
                           by = "ProteinName"]
        merge(unique(input_loop[, .(ProteinName, Run, Channel, CenteredAbundance = 0)]),
              means, by = "ProteinName")
    } else {
        input_loop[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
        prot_has_unique = input_loop[(IsUnique), unique(ProteinName)]
        initial_input = input_loop[ProteinName %in% prot_has_unique]
        initial_weights = unique(initial_input[ProteinName %in% prot_has_unique,
                                               .(ProteinName, PSM)])
        initial_weights[, Weight := 1 / uniqueN(ProteinName), by = "PSM"]
        initial_unique = summarizeProteinsClusterSingleRun(initial_input,
                                                           initial_weights,
                                                           norm, norm_parameter,
                                                           use_shared = FALSE)
        means = input_loop[!(ProteinName %in% prot_has_unique),
                           .(Abundance = mean(log2IntensityNormalized,
                                              na.rm = TRUE)),
                           by = "ProteinName"]
        initial_no_unique = merge(unique(input_loop[, .(ProteinName, Run,
                                                        Channel,
                                                        CenteredAbundance = 0)]),
                                  means, by = "ProteinName")
        rbind(initial_unique, initial_no_unique)
    }

}


#' Utility function to get weights from current iterations including 0s (removed features)
#' @param weights table of weights (columns ProteinName, PSM, Weight)
#' @inheritParams getWeightedSummarySingleRun
#' @keywords internal
getCurrentWeights = function(weights, peptide_protein_dt) {
    weights = merge(weights, peptide_protein_dt,
                    all.x = TRUE, all.y = TRUE,
                    by = c("ProteinName", "PSM"))
    weights[, Weight := ifelse(is.na(Weight), 0, Weight)]
    weights[order(PSM, ProteinName)][, Weight]
}
