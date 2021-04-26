#' Plot multiple protein-level summaries
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
            geom_boxplot(aes(x = Channel, y = log2IntensityNormalized),
                         data = feature_data, inherit.aes = FALSE) +
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
        theme(axis.text.x = element_text(angle = 270))
    plot
}


#' Compare distribution of log2-intensities per channel with and without shared peptides
#'
#' @param input data.table
#' @param plot_profiles if TRUE, feature profiles will be plotted
#'
#' @return ggplot
#'
#' @export
#'
plotDistributionComparison = function(input, plot_profiles = FALSE) {
    full_distribution = input[, .(ProteinName, Run, PSM, Channel, log2IntensityNormalized,
                                  IsUnique, Type = "All peptides")]
    unique_distribution = input[(IsUnique), .(ProteinName, Run, PSM, Channel, log2IntensityNormalized,
                                              IsUnique, Type = "Unique peptides")]
    distributions = rbind(full_distribution, unique_distribution)
    distributions$IsUnique = factor(as.character(distributions$IsUnique),
                                    levels = c("TRUE", "FALSE"),
                                    ordered = TRUE)
    plot = ggplot(distributions, aes(x = Channel, y = log2IntensityNormalized,
                                     fill = Type)) +
        geom_boxplot()
    if (plot_profiles) {
        plot = plot +
            geom_line(aes(x = Channel, y = log2IntensityNormalized,
                          linetype = IsUnique, group = PSM),
                      color = "grey", size = 0.8, alpha = 0.6)
    }
    plot +
        facet_wrap(Run ~ ProteinName) +
        theme_bw()
}


#' Plot feature profiles
#'
#' @param input data.table
#'
#' @return ggplot
#'
#' @export
#'
plotProfiles = function(input) {
    input$IsUnique = factor(as.character(input$IsUnique),
                            levels = c("TRUE", "FALSE"), ordered = TRUE)
    ggplot(input, aes(x = Channel, y = log2IntensityNormalized,
                      group = PSM, linetype = IsUnique)) +
        geom_point() +
        geom_line() +
        facet_wrap(Run ~ ProteinName) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 270))
}
