feature_data = readRDS("../SharedPeptideTMTReproduction/brd_triple_example.RDS")
norm = "Huber"
norm_parameter = 0.1
weights_mode = "contributions"
use_control = TRUE
control_pattern = NULL
tolerance = 1e-3
max_iter = 10

brd_summ = getWeightedProteinSummary(feature_data, norm, norm_parameter,
                                     max_iter = 30, tolerance = 0.03)
brd_summ$ConvergenceSummary

brd_summ$Weights

summ = brd_summ$ProteinLevelData
ch_order = c("126C", "127N", "127C", "128N", "128C",
             "129N", "129C", "130N", "130C", "131N")
summ[, Channel := factor(Channel, levels = ch_order, ordered = TRUE)]
feats = brd_summ$FeatureLevelData
feats[, Channel := factor(Channel, levels = ch_order, ordered = TRUE)]

unique_summary = summarizeProteinsClusterSingleRun(feature_data,
                                                   unique(feature_data[, .(ProteinName, PSM, Weight = 1)]),
                                                   "Huber", 0.1, FALSE)

all_summary = data.table::rbindlist(lapply(
    split(feature_data, feature_data[, ProteinName]),
    function(x) {
        summarizeProteinsClusterSingleRun(x,
                                          unique(x[, .(ProteinName, PSM, Weight = 1)]),
                                          "Huber", 0.1, FALSE)
    }
))
all_summary[, Channel := factor(Channel, levels = ch_order, ordered = TRUE)]



library(ggplot2)
ggplot(summ, aes(x = Channel, y = Abundance, group = ProteinName,
                 color = ProteinName)) +
    geom_line(aes(x = Channel, y = log2IntensityNormalized,
                  group = PSM),
              data = feats, inherit.aes = FALSE) +
    geom_line(aes(x = Channel, y = Abundance, group = ProteinName,
                  color = ProteinName),
              data = all_summary, color = "black", size = 1.3) +
    geom_line(size = 1.5) +
    facet_grid(~ProteinName) +
    theme_bw() +
    theme(legend.position = "bottom")


brd_summ$ProteinLevelData




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

library(ggplot2)
ggplot(summarized_output, aes(x = Channel, y = Abundance, group = ProteinName, color = ProteinName)) +
    geom_line() +
    theme_bw()



summarizeProteinsClusterSingleRun = function(input, peptide_protein_dt,
                                             norm, norm_parameter, weights_mode,
                                             use_control, control_pattern,
                                             tolerance, max_iter) {
    alphas_list = vector("list", max_iter)
    input_loop = data.table::copy(input)
    input_loop = unique(input_loop[, .(ProteinName, PSM, Channel, log2IntensityNormalized)])
    num_alpha_params = data.table::uniqueN(peptide_protein_dt)
    current_alphas = rep(1, num_alpha_params)
    previous_alphas = rep(0, num_alpha_params)

    initial_summary_input = data.table::copy(input_loop)
    weights_df = unique(initial_summary_input[, list(ProteinName, PSM, Weight = 1)])
    initial_summary = summarizeProteinsCluster(initial_summary_input, weights_df, norm,
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
                                  use_control, control_pattern)
        alphas_df_round = alphas_df[Weight > 0.01]
        alphas_df_round[, Weight := Weight / sum(Weight), by = "PSM"]
        alphas_list[[iter]] = alphas_df_round
        input_loop = merge(input_loop,
                           unique(alphas_df_round[, .(PSM, ProteinName)]),
                           by = c("ProteinName", "PSM"))

        new_abundances = summarizeProteinsCluster(input_loop, alphas_df_round,
                                                  norm, norm_parameter, TRUE)

        input_loop = merge(input_loop[, list(ProteinName, PSM, Channel,
                                             log2IntensityNormalized)],
                           new_abundances, by = c("ProteinName", "Channel"))
        alphas_df = merge(alphas_df, peptide_protein_dt,
                          all.x = TRUE, all.y = TRUE,
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
    summarized_output = new_abundances
    summarized_output[, Abundance := Abundance + global_intercept]
    summarized_output = summarized_output[, !(colnames(summarized_output) %in% c("PSM", "log2IntensityNormalized", "global_intercept")), with = FALSE]
    summarized_output = unique(summarized_output)
    alpha_diffs = alpha_diffs[alpha_diffs > 0]
    alphas_list = alphas_list[sapply(alphas_list, function(x) !is.null(x))]
    list(summary = summarized_output,
         alpha_history = alphas_list,
         convergence_history = alpha_diffs)
}

