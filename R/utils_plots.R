#' Plot protein-level summaries for a selected cluster of proteins
#'
#' @param weighted_summary Output of the getWeightedProteinSummary function
#' @param cluster ID of a cluster of proteins summarized in the `weighted_summary` object
#' @param channel_order optional character vector of labels of Runs or Channels which
#' will be used to sort the x-axis
#'
#' @return ggplot2 object
#'
#' @export
#'
plotSummarizedProteins = function(weighted_summary, cluster, channel_order = NULL) {
    feature_plot_input = featureData(weighted_summary)[Cluster == cluster]
    feature_plot_input[, IsUnique := data.table::uniqueN(ProteinName) == 1,
                       by = "PSM"]
    feature_plot_input[, Peptide := factor(ifelse(IsUnique, "unique", "shared"),
                                           levels = c("unique", "shared"),
                                           ordered = TRUE)]

    if (!is.null(channel_order)) {
        feature_plot_input[, Channel := factor(Channel, levels = channel_order,
                                               ordered = TRUE)]
        protein_plot_input[, Channel := factor(Channel, levels = channel_order,
                                               ordered = TRUE)]
    }

    protein_plot_input = proteinData(weighted_summary)[Protein %in% unique(feature_plot_input$ProteinName)]
    setnames(protein_plot_input, "Protein", "ProteinName")
    if (weighted_summary@ExperimentType == "LF") {
        x_axis = "Run"
        annot = "run"
    } else {
        x_axis = "Channel"
        annot = "channel"
    }

    plot = ggplot() +
        geom_line(aes_string(x = x_axis, y = "log2IntensityNormalized",
                             group = "PSM", linetype = "Peptide"),
                  data = feature_plot_input, alpha = 0.5, size = 1.2) +
        geom_line(aes_string(x = x_axis, y = "Abundance",
                             color = "ProteinName", group = "ProteinName"),
                  data = protein_plot_input, size = 2) +
        scale_linetype_discrete(name = "peptide") +
        scale_color_discrete(palette = "viridis") +
        xlab(annot) +
        ylab("log-intensity") +
        guides(color = "none") +
        theme_bw() +
        theme(legend.position = "bottom",
              axis.text.x = element_text(angle = 270),
              legend.direction = "horizontal")
    if (weighted_summary@ExperimentType == "LF") {
        plot = plot +
            facet_grid( ~ ProteinName)
    } else {
        plot = plot +
            facet_grid(Run ~ ProteinName)
    }
    plot
}

#' Plot multiple protein-level summaries for the same set of proteins
#'
#' @param ... data.tables with summaries
#' @param channel_order optional order for x-axis (Channel column)
#' @param feature_data optional data.table for plotting feature-level distribution
#' in each channel
#'
#' @return ggplot
#' @import ggplot2
#'
#' @export
#'
plotSummaryComparison = function(..., channel_order = NULL, feature_data = NULL) {
    df_list = list(...)
    summaries_df = rbindlist(df_list, use.names = TRUE, fill = TRUE)
    if (!is.null(channel_order)) {
        summaries_df$Channel = factor(summaries_df$Channel,
                                      levels = channel_order,
                                      ordered = TRUE)
    }
    plot = ggplot(summaries_df, aes(x = Channel, y = Abundance, group = Method,
                                    color = Method))
    if (!is.null(feature_data)) {
        plot = plot +
            # geom_boxplot(aes(x = Channel, y = log2IntensityNormalized),
            #              data = feature_data, inherit.aes = FALSE) +
            geom_line(aes(x = Channel, y = log2IntensityNormalized,
                          group = PSM, linetype = IsUnique),
                      data = feature_data, inherit.aes = FALSE,
                      color = "grey", alpha = 0.5, size = 0.8)
    }
    plot = plot +
        geom_point(size = 1.2) +
        geom_line(size = 1.2) +
        facet_wrap(Run~ProteinName) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 270),
              legend.position = "bottom")
    plot
}


#' Plot observed PSM-level profiles and compare them to profiles predicted by weighted summarization model
#'
#' @inheritParams plotSummarizedProteins
#'
#' @return ggplot2
#'
#' @export
#'
plotFittedProfiles = function(weighted_summary, cluster, channel_order = NULL) {
    fitted_profiles = fittedProfiles(weighted_summary)[Cluster == cluster]
    fitted_profiles = melt(fitted_profiles,
                           measure.vars = c("log2IntensityNormalized",
                                            "Predicted"),
                           variable.factor = FALSE)
    fitted_profiles[, Profile := ifelse(variable == "Predicted",
                                        "fitted", "observed")]
    fitted_profiles[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
    fitted_profiles[, Peptide := ifelse(IsUnique, "unique", "shared")]
    fitted_profiles[, Peptide := factor(Peptide,
                                        levels = c("unique", "shared"),
                                        ordered = TRUE)]
    if (weighted_summary@ExperimentType == "TMT") {
        if (!is.null(channel_order)) {
            fitted_profiles[, Channel := factor(Channel, levels = channel_order,
                                                ordered = TRUE)]
        }
    } else {
        if (!is.null(channel_order)) {
            fitted_profiles[, Run := factor(Run, levels = channel_order,
                                            ordered = TRUE)]
        }
    }
    if (weighted_summary@ExperimentType == "LF") {
        x_axis = "Run"
        annot = "run"
    } else {
        x_axis = "Channel"
        annot = "channel"
    }

    fitted_profiles[, grouping := paste(variable, PSM)]
    plot = ggplot(fitted_profiles, aes_string(x = x_axis, y = "value",
                                              group = "grouping",
                                              color = "Profile",
                                              linetype = "Peptide")) +
        geom_line(linewidth = 1.2) +
        scale_color_discrete(name = "fitted", palette = "viridis") +
        scale_linetype_discrete(name = "peptide") +
        xlab("annot") +
        ylab("log-intensity") +
        theme_bw() +
        theme(legend.position = "bottom",
              axis.text.x = element_text(angle = 270))
    if (weighted_summary@ExperimentType == "LF") {
        plot = plot +
            facet_grid( ~ ProteinName)
    } else {
        plot = plot +
            facet_grid(Run ~ ProteinName)
    }
    plot
}
