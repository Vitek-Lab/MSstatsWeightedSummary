#' Wrapper around `igraph::graph_from_data_frame`
#'
#' @param quantification_data MS data, preferably in `MSstats` or `MSstatsTMT` format.
#' @param protein_column name of a column with protein names.
#' @param peptide_column name of a column with peptide sequences.
#' @param by_run applies only to TMT data. If TRUE, clusters will be computed for each TMT run separately.
#'
#' @importFrom igraph graph_from_data_frame
#'
#' @export
#'
createPeptideProteinGraph = function(quantification_data,
                                     protein_column = "ProteinName",
                                     peptide_column = "PeptideSequence",
                                     by_run = FALSE) {
    if (by_run) {
        lapply(split(quantification_data, quantification_data[["Run"]]),
               function(single_run) {
                   igraph::graph_from_data_frame(
                       unique(single_run[, c(protein_column, peptide_column),
                                         with = FALSE]),
                       directed = FALSE)
               })
    } else {
        igraph::graph_from_data_frame(
            unique(quantification_data[, c(protein_column, peptide_column),
                                       with = FALSE]),
            directed = FALSE)
    }

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
    Run = ProteinName = NULL

    if (inherits(peptide_protein_graph, "list")) {
        membership = data.table::rbindlist(lapply(names(peptide_protein_graph), function(run_name) {
                all_unique_proteins = unique(quantification_data[Run == run_name, ProteinName])
                graph_decomposed = igraph::decompose.graph(peptide_protein_graph[[run_name]])
                data.table::rbindlist(lapply(
                    1:length(graph_decomposed),
                    function(cluster_id) {
                        list(Cluster = paste(run_name, cluster_id, sep = "__"),
                             Run = run_name,
                             ProteinName = intersect(names(igraph::V(graph_decomposed[[cluster_id]])),
                                                     all_unique_proteins))
                    }))
            }))
        merge(quantification_data, membership, by = c("ProteinName", "Run"),
              all.x = TRUE, sort = FALSE)
    } else {
        all_unique_proteins = unique(quantification_data[[protein_column]])
        graph_decomposed = igraph::decompose.graph(peptide_protein_graph)
        membership = data.table::rbindlist(lapply(
            1:length(graph_decomposed),
            function(cluster_id) {
                data.table::data.table(Cluster = cluster_id,
                                       ProteinName = intersect(names(igraph::V(graph_decomposed[[cluster_id]])),
                                                               all_unique_proteins))
            }))
        merge(quantification_data, membership, by = "ProteinName",
              all.x = TRUE, sort = FALSE)
    }
}


#' Calculate statistics that describe clusters of proteins and peptides.
#'
#' @param quantification_data MS data, preferably in `MSstats` or `MSstatsTMT` format.
#' @param merge if `TRUE`, calculated statistics will be merged into original data.
#'
#' @export
#'
getClusterStatistics = function(quantification_data, merge = FALSE) {
    `:=` = NumProteins = NumPeptides = ProteinName = PeptideSequence = EachHasUnique = NULL
    NumProteinsPerPeptide = TotalSize = NumPeptidesPerProtein = IsUnique = HasUnique = AnyHasUnique = NULL

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
                           by = c("Cluster", "ProteinName", "PeptideSequence"),
                           sort = FALSE)
    }
    statistics
}
