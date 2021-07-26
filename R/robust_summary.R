#' Get robust summary with shared peptides
#'
#' @param input data.table in MSstatsTMT format
#' @param method_label name for a summarization method
#' @param design_type "short" or "long" (not implemented yet)
#' @param norm "p_norm" or "Huber"
#' @param norm_parameter p for norm=="p_norm", M for norm=="Huber"
#' @param use_weights if TRUE, PSM-Protein weights will be used
#' @param weights_method if "separate", weights will be estimated for each PSM
#' separately, if "all", weights will be estimated jointly
#' @param adaptive_huber if TRUE, M parameter for Huber norm will be updated
#' in each iteration of the procedure
#' @param initial_summary_method summarization method to start iterative procedure
#' - "msstats" or "same"
#' @param use_discordant if TRUE, discordant profiles will be used
#' @param use_shared_initial if TRUE, shared peptides will be used in summary
#' that initializes iterative procedure
#' @param tolerance tolerance to indicate weights convergence
#' @param max_iter maximum number of iteration of the procedure
#'
#' @return list of final summary, history of weights and difference between weights
#' from consecutive iterations
#'
#' @export
#'
getRobustSummary = function(input,
                            method_label= "",
                            norm = "p_norm", norm_parameter = 1,
                            weights_method = "separate",
                            adaptive_huber = TRUE,
                            initial_summary_method = "msstats",
                            use_discordant = TRUE,
                            equalize_protein_features = FALSE,
                            tolerance = 1e-3,
                            max_iter = 10) {
    if (!is.element("Weight", colnames(input))) {
        input[, Weight := 1 / uniqueN(ProteinName), by = "PSM"]
    }
    if (any(!is.finite(input$log2IntensityNormalized) | is.na(input$log2IntensityNormalized))) {
        input[,
              log2IntensityNormalized := ifelse(is.na(log2IntensityNormalized) |
                                                    !is.finite(log2IntensityNormalized),
                                                median(log2IntensityNormalized, na.rm = TRUE),
                                                log2IntensityNormalized), by = c("ProteinName", "Run", "Channel")]
    }
    input_by_run = split(input, input$Run)
    output_by_run = lapply(input_by_run, function(x) {
        getRobustSummarySingleRun(x, method_label,
                                  norm, norm_parameter,
                                  weights_method, adaptive_huber,
                                  initial_summary_method, use_discordant,
                                  equalize_protein_features,
                                  tolerance, max_iter)
    })
    summarized_output = rbindlist(lapply(output_by_run, function(x) x$summary))
    alphas_list = lapply(output_by_run, function(x) x$alpha_history)
    alpha_diffs = lapply(output_by_run, function(x) x$convergence_history)
    list(summary = summarized_output,
         alpha_history = alphas_list,
         convergence_history = alpha_diffs)
}


#' Robust summary for a single run
#' @keywords internal
getRobustSummarySingleRun = function(input,
                                     method_label,
                                     norm = "p_norm", norm_parameter = 1,
                                     weights_method = "all",
                                     adaptive_huber = TRUE,
                                     initial_summary_method = "msstats",
                                     use_discordant = TRUE,
                                     equalize_protein_features = FALSE,
                                     tolerance = 1e-3,
                                     max_iter = 10) {
    annotation = unique(input[, list(Run, Mixture, TechRepMixture,
                                     Channel, Condition, BioReplicate)])
    alphas_list = vector("list", max_iter)
    input_loop = copy(input)
    input_loop = unique(input_loop[, .(ProteinName, PSM, Channel, log2IntensityNormalized, Weight)])
    ppsm = unique(input[, .(ProteinName, PSM)])
    num_alpha_params = uniqueN(ppsm)
    current_alphas = rep(1, num_alpha_params)
    previous_alphas = rep(0, num_alpha_params)
    initial_summary = getInitialSummary(input, method = initial_summary_method,
                                        use_discordant, use_shared_initial,
                                        norm, norm_parameter)
    input_loop = merge(input_loop,
                       initial_summary[, list(ProteinName, Channel, Abundance)],
                       by = c("ProteinName", "Channel"))
    iter = 1
    alpha_diffs = vector("numeric", max_iter)
    while (sum(abs(current_alphas - previous_alphas)) > tolerance) {
        previous_alphas = current_alphas
        alphas_df = getAllWeights(input_loop, weights_method, norm, norm_parameter,
                                  equalize_protein_features)
        alphas_df_round = alphas_df[Weight > 0.01]
        alphas_df_round[, Weight := Weight / sum(Weight), by = "PSM"]
        alphas_list[[iter]] = alphas_df_round

        input_abundances = merge(input_loop[, list(ProteinName, PSM, Channel,
                                                   log2IntensityNormalized)],
                                 alphas_df_round, by = c("ProteinName", "PSM"))
        new_abundances = getOptimSummary(input_abundances, "short", TRUE,
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
    summarized_output = summarized_output[, !(colnames(summarized_output) %in% c("PSM", "log2IntensityNormalized")), with = FALSE]
    summarized_output = unique(summarized_output)
    summarized_output$Method = method_label
    alpha_diffs = alpha_diffs[alpha_diffs > 0]
    alphas_list = alphas_list[sapply(alphas_list, function(x) !is.null(x))]
    list(summary = summarized_output,
         alpha_history = alphas_list,
         convergence_history = alpha_diffs)
}

