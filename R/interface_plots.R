#' Plot protein-level summaries for a selected cluster of proteins
#'
#' @param weighted_summary Output of the getWeightedProteinSummary function
#' @param cluster ID of a cluster of proteins summarized in the `weighted_summary` object
#' @param channel_order optional character vector of labels of Runs or Channels which
#' will be used to sort the x-axis
#'
#' @return ggplot2 object
#' @import ggplot2
#'
#' @export
#'
plotSummarizedProteins = function(weighted_summary, cluster, channel_order = NULL) {
    `:=` = Cluster = IsUnique = ProteinName = Peptide = Protein = Channel = NULL

    feature_plot_input = featureData(weighted_summary)[Cluster == cluster]
    feature_plot_input[, IsUnique := data.table::uniqueN(ProteinName) == 1,
                       by = "PSM"]
    feature_plot_input[, Peptide := factor(ifelse(IsUnique, "unique", "shared"),
                                           levels = c("unique", "shared"),
                                           ordered = TRUE)]

    protein_plot_input = proteinData(weighted_summary)[Protein %in% unique(feature_plot_input[["ProteinName"]])]
    data.table::setnames(protein_plot_input, "Protein", "ProteinName")

    if (!is.null(channel_order)) {
        feature_plot_input[, Channel := factor(Channel, levels = channel_order,
                                               ordered = TRUE)]
        protein_plot_input[, Channel := factor(Channel, levels = channel_order,
                                               ordered = TRUE)]
    }

    if (weighted_summary@ExperimentType == "LF") {
        x_axis = "Run"
        annot = "run"
    } else {
        x_axis = "Channel"
        annot = "channel"
    }

    plot = ggplot() +
        geom_line(aes(x = .data[[x_axis]], y = .data[["log2IntensityNormalized"]],
                             group = .data[["PSM"]], linetype = .data[["Peptide"]]),
                  data = feature_plot_input, alpha = 0.5, linewidth = 1.2) +
        geom_line(aes(x = .data[[x_axis]], y = .data[["Abundance"]],
                             color = .data[["ProteinName"]], group = .data[["ProteinName"]]),
                  data = protein_plot_input, linewidth = 2) +
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


#' Plot observed PSM-level profiles and compare them to profiles predicted by weighted summarization model
#'
#' @inheritParams plotSummarizedProteins
#'
#' @return ggplot2
#'
#' @export
#'
plotFittedProfiles = function(weighted_summary, cluster, channel_order = NULL) {
    `:=` = Cluster = Profile = IsUnique = Peptide = Channel = Run = PSM = variable = ProteinName = NULL

    fitted_profiles = fittedProfiles(weighted_summary)[Cluster == cluster]
    fitted_profiles = data.table::melt(fitted_profiles,
                           measure.vars = c("log2IntensityNormalized",
                                            "Predicted"),
                           variable.factor = FALSE)
    fitted_profiles[, Profile := ifelse(variable == "Predicted",
                                        "fitted", "observed")]
    fitted_profiles[, IsUnique := data.table::uniqueN(ProteinName) == 1, by = "PSM"]
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
    plot = ggplot(fitted_profiles, aes(x = .data[[x_axis]], y = .data[["value"]],
                                              group = .data[["grouping"]],
                                              color = .data[["Profile"]],
                                              linetype = .data[["Peptide"]])) +
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
