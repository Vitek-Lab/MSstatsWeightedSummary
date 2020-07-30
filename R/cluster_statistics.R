#' Wrapper around `igraph::graph_from_data_frame`
#'
#' @param quantification_data MS data, preferably in `MSstats` or `MSstatsTMT` format.
#' @param protein_column name of a column with protein names.
#' @param peptide_column name of a column with peptide sequences.
#'
#' @importFrom igraph graph_from_data_frame
#'
#' @export
#'
createPeptideProteinGraph = function(quantification_data,
                                   protein_column = "ProteinName",
                                   peptide_column = "PeptideSequence") {
    igraph::graph_from_data_frame(
        unique(quantification_data[, c(protein_column, peptide_column),
                                   with = FALSE]),
        directed = FALSE)
}

#' Add information about connected subcomponents of the peptide-protein graph
#' to quantitative data.
#'
#' @param quantification_data MS data, preferably in `MSstats` or `MSstatsTMT` format.
#' @param peptide_protein_graph graph created by the `createPeptideProteinGraph` function.
#' @param protein_column name of a column with protein names.
#'
#' @importFrom igraph decompose.graph V
#' @importFrom data.table data.table
#'
#' @export
#'
addClusterMembership = function(quantification_data, peptide_protein_graph,
                                protein_column = "ProteinName") {
    graph_decomposed = igraph::decompose.graph(peptide_protein_graph)
    all_unique_proteins = unique(quantification_data[[protein_column]])

    membership = data.table::rbindlist(lapply(
        1:length(graph_decomposed),
        function(cluster_id) {
            data.table::data.table(Cluster = cluster_id,
                                   ProteinName = intersect(names(igraph::V(graph_decomposed[[cluster_id]])),
                                                           all_unique_proteins))
        }))
    merge(quantification_data, membership, by = "ProteinName", all.x = TRUE)
}


#' Calculate statistics that describe clusters of proteins and peptides.
#'
#' @param quantification_data MS data, preferably in `MSstats` or `MSstatsTMT` format.
#' @param merge if `TRUE`, calculated statistics will be merged into original data.
#'
#' @export
#'
getClusterStatistics = function(quantification_data, merge = FALSE) {
    ProteinName = PeptideSequence = NULL

    statistics = quantification_data[, list(ProteinName, PeptideSequence,
                                            NumProteins = data.table::uniqueN(ProteinName),
                                            NumPeptides = data.table::uniqueN(PeptideSequence)),
                                     by = "Cluster"]
    statistics = unique(statistics)
    statistics[, TotalSize := NumProteins + NumPeptides]
    statistics[, NumProteinsPerPeptide := data.table::uniqueN(ProteinName),
               by = "PeptideSequence"]
    statistics[, NumPeptidesPerProtein := data.table::uniqueN(PeptideSequence),
               by = "ProteinName"]
    statistics[, IsUnique := NumProteinsPerPeptide == 1L]
    statistics[, HasUnique := any(IsUnique), by = "ProteinName"]
    statistics[, AnyHasUnique := any(HasUnique), by = "Cluster"]
    statistics[, EachHasUnique := all(HasUnique), by = "Cluster"]


    if (merge) {
        statistics = merge(quantification_data, statistics,
                           by = c("Cluster", "ProteinName", "PeptideSequence"))
    }
    statistics
}

#' Plot statistics about clusters of proteins
#'
#' @param cluster_statistics output of the `getClusterStatistics` function.
#' @param statistic one of the column names of the `cluster_statistics` parameter.
#' @param only_clusters_with_shared if `TRUE`, statistics will be presented only
#' for cluster with more than one protein.
#'
#' @import ggplot2
#'
#' @export
#'
plotClusterStats = function(cluster_statistics, statistic, only_clusters_with_shared = FALSE) {
    if (only_clusters_with_shared) {
        shared_filter = cluster_statistics$NumProteins > 1
    } else {
        shared_filter = rep(TRUE, nrow(cluster_statistics))
    }

    if (statistic %in% c("NumProteins", "NumPeptides", "TotalSize",
                         "MedProteinsPerPeptide", "MedPeptidesPerProtein",
                         "AnyHasUnique", "EachHasUnique")) {
        data_to_plot = unique(cluster_statistics[shared_filter, c("Cluster", statistic), with = FALSE])
        y_label = "number of clusters"
    } else if (statistic %in% c("NumPeptidesPerProtein", "HasUnique")) {
        data_to_plot = unique(cluster_statistics[shared_filter, c("ProteinName", statistic), with = FALSE])
        y_label = "number of proteins"
    } else if (statistic %in% c("NumProteinsPerPeptide", "IsUnique")) {
        data_to_plot = unique(cluster_statistics[shared_filter, c("PeptideSequence", statistic), with = FALSE])
        y_label = "number of peptides"
    }

    if (is.logical(data_to_plot[[2]])) {
        geom = geom_bar()
    } else {
        geom = geom_histogram(binwidth = 1)
    }

    ggplot(data_to_plot, aes_string(x = statistic)) +
        geom +
        # scale_y_log10() +
        ylab(y_label) +
        theme_bw()
}
