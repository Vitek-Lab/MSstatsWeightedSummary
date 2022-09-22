#' Get robust summary with shared peptides
#'
#' @param input data.table in MSstatsTMT format
#' @param method_label name for a summarization method
#' @param norm "p_norm" or "Huber"
#' @param norm_parameter p for norm=="p_norm", M for norm=="Huber"
#' @param weights_mode "contributions" for "sum to one" condition,
#' "probabilities" for only "non-negative" condition.
#' @param tolerance tolerance to indicate weights convergence
#' @param max_iter maximum number of iteration of the procedure
#'
#' @return list of final summary, history of weights and difference between weights
#' from consecutive iterations
#'
#' @export
#'
getRobustSummary = function(input,
                            norm = "p_norm", norm_parameter = 1,
                            weights_mode = "contributions",
                            use_control = TRUE, control_pattern = NULL,
                            tolerance = 1e-3, max_iter = 10) {
    input = data.table::as.data.table(input)
    input_by_run = split(input, input[, Run])

    output_by_run = lapply(input_by_run, function(x) {
        getRobustSummarySingleRun(x, norm, norm_parameter,
                                  weights_mode, use_control,
                                  control_pattern, tolerance, max_iter)
    })

    summarized_output = data.table::rbindlist(lapply(output_by_run,
                                                     function(x) x[["summary"]]))
    alphas_list = lapply(output_by_run, function(x) x[["alpha_history"]])
    alpha_diffs = lapply(output_by_run, function(x) x[["convergence_history"]])

    list(summary = summarized_output,
         alpha_history = alphas_list,
         convergence_history = alpha_diffs,
         tolerance = tolerance)
}


#' Robust summary for a single run
#' @keywords internal
getRobustSummarySingleRun = function(input,
                                     norm = "p_norm", norm_parameter = 1,
                                     weights_mode = "contributions",
                                     use_control = TRUE, control_pattern = NULL,
                                     tolerance = 1e-3, max_iter = 10) {

    annotation = unique(input[, list(Run, Mixture, TechRepMixture,
                                     Channel, Condition, BioReplicate)])
    alphas_list = vector("list", max_iter)
    input_loop = copy(input)
    input_loop = unique(input_loop[, .(ProteinName, PSM, Channel, log2IntensityNormalized)])
    ppsm = unique(input[, .(ProteinName, PSM)])
    num_alpha_params = uniqueN(ppsm)
    current_alphas = rep(1, num_alpha_params)
    previous_alphas = rep(0, num_alpha_params)

    initial_summary_input = data.table::copy(input_loop)
    weights_df = unique(initial_summary_input[, list(ProteinName, PSM, Weight = 1)])
    initial_summary = getWeightedSummary(initial_summary_input, weights_df, norm,
                                         norm_parameter, FALSE)

    input_loop = merge(input_loop,
                       initial_summary[, list(ProteinName, Channel, Abundance)],
                       by = c("ProteinName", "Channel"))
    iter = 1
    alpha_diffs = vector("numeric", max_iter)
    while (sum(abs(current_alphas - previous_alphas)) > tolerance) {
        previous_alphas = current_alphas
        alphas_df = getAllWeights(input_loop, norm,
                                  norm_parameter, weights_mode,
                                  use_control = TRUE, control_pattern = NULL)
        alphas_df_round = alphas_df[Weight > 0.01]
        alphas_df_round[, Weight := Weight / sum(Weight), by = "PSM"]
        alphas_list[[iter]] = alphas_df_round
        input_loop = merge(input_loop,
                           unique(alphas_df_round[, .(PSM, ProteinName)]),
                           by = c("ProteinName", "PSM"))

        new_abundances = getWeightedSummary(input_loop, alphas_df_round,
                                            norm, norm_parameter, TRUE,
                                            adaptive_huber)
        if (norm == "Huber" & adaptive_huber) {
            norm_parameter = attr(new_abundances, "M")
        }
        input_loop = merge(input_loop[, list(ProteinName, PSM, Channel,
                                             log2IntensityNormalized)],
                           new_abundances, by = c("ProteinName", "Channel"))
        alphas_df = merge(alphas_df, ppsm, all.x = TRUE, all.y = TRUE,
                          by = c("ProteinName", "PSM"))
        alphas_df[, Weight := ifelse(is.na(Weight), 0, Weight)]
        current_alphas = alphas_df$Weight

        if (iter >= max_iter) {
            current_alphas = previous_alphas
        } else {
            alpha_diffs[iter] = sum(abs(current_alphas - previous_alphas))
            iter = iter + 1
        }
    }
    summarized_output = merge(input_loop, annotation,
                              by = "Channel")
    summarized_output[, Abundance := Abundance + global_intercept]
    summarized_output = summarized_output[, !(colnames(summarized_output) %in% c("PSM", "log2IntensityNormalized", "global_intercept")), with = FALSE]
    summarized_output = unique(summarized_output)
    summarized_output$Method = method_label
    alpha_diffs = alpha_diffs[alpha_diffs > 0]
    alphas_list = alphas_list[sapply(alphas_list, function(x) !is.null(x))]
    list(summary = summarized_output,
         alpha_history = alphas_list,
         convergence_history = alpha_diffs)
}

