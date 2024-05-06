#' @export
processIsoforms = function(quantification_data, remove_single_shared = TRUE,
                           merge_identical = TRUE, remove_subsets = FALSE,
                           subset_treatment = "remove_all") {
    quantification_data = getUniquenessInfo(quantification_data)
    unique_only = quantification_data[(UniqueOnly)]
    quantification_data = quantification_data[!(UniqueOnly)]

    if (remove_single_shared) {
        quantification_data = quantification_data[!(NumPeptidesPerProtein == 1 & !(HasUnique))]
        quantification_data = getUniquenessInfo(quantification_data)
        unique_only = rbind(unique_only, quantification_data[(UniqueOnly)])
        quantification_data = quantification_data[!(UniqueOnly)]
    }

    if (merge_identical) {
        pp_graph = createPeptideProteinGraph(quantification_data)
        quantification_data[, Cluster := NULL]
        quantification_data = addClusterMembership(quantification_data,
                                                   pp_graph)
        quantification_data[, Cluster := Cluster + max(unique_only$Cluster)]
        data_by_cluster = split(quantification_data, quantification_data[["Cluster"]])
        processed_clusters = lapply(data_by_cluster, mergeIdenticalProteins)
        processed_clusters = data.table::rbindlist(processed_clusters)
        setnames(quantification_data, "ProteinName", "ProteinNameOriginal")
        quantification_data = merge(quantification_data, processed_clusters,
                                    by = "ProteinNameOriginal")
        quantification_data = unique(quantification_data[, colnames(quantification_data) != "ProteinNameOriginal", with = FALSE])
        quantification_data = getUniquenessInfo(quantification_data)
    }

    if (remove_subsets) {
        if (subset_treatment == "remove_all") {
            quantification_data = quantification_data[(HasUnique)]
            # Below can be removed?
            quantification_data = getUniquenessInfo(quantification_data)
            unique_only = rbind(unique_only, quantification_data[(UniqueOnly)],
                                fill = TRUE)
            quantification_data = quantification_data[!(UniqueOnly)]
        } else {
            subset_clusters = quantification_data[!(HasUnique)]
            subset_counts = subset_clusters[, .(NumFeatures = data.table::uniqueN(PSM)),
                                            by = c("ProteinName", "Cluster")]
            subset_counts[, Rank := rank(-NumFeatures), by = "Cluster"]
            subset_counts[, list(ProteinName = paste(sort(unique(ProteinName)),
                                                     sep = ";", collapse = ";"),
                                 OriginalProteinName = ProteinName[is.max(Rank)]), by = "Cluster"]
            # subset_counts = subset_counts[Rank == 1]
            # quantification_data = quantification_data[(HasUnique) | ProteinName %in% subset_counts[, ProteinName]]
            # can be done with merge?
            quantification_data[, ProteinName := NULL]
            quantification_data = merge(quantificatioN_data, subset_counts, )
            quantification_data = getUniquenessInfo(quantification_data)
        }
    }

    output = rbind(unique_only, quantification_data, fill = TRUE)
    output = output[, .(ProteinName, PeptideSequence, Charge, PSM,
                        Run, Channel, log2Intensity, BioReplicate, Condition,
                        Mixture, TechRepMixture)]
    pp_g = createPeptideProteinGraph(output)
    output = addClusterMembership(output, pp_g)
    output
}

#' @keywords internal
getUniquenessInfo = function(quantification_data) {
    quantification_data[, IsUnique := data.table::uniqueN(ProteinName) == 1,
                        by = "PSM"]
    quantification_data[, HasUnique := any(IsUnique), by = "ProteinName"]
    quantification_data[, NumPeptidesPerProtein := data.table::uniqueN(PSM),
                        by = "ProteinName"]
    quantification_data[, UniqueOnly := all(IsUnique), by = "ProteinName"]
    quantification_data
}

#' @keywords internal
mergeIdenticalProteins = function(cluster_data) {
    pp_m = as.matrix(data.table::dcast(unique(cluster_data[, .(ProteinName, PSM, Present = 1)]),
                                     PSM ~ ProteinName, value.var = "Present",
                                     fill = 0)[, -1])
    pp_m_cor = getIdenticalIndicator(pp_m)
    iso_graph = igraph::graph_from_adjacency_matrix(abs(pp_m_cor - 1) < 1e-10)
    graph_decomposed = igraph::decompose.graph(iso_graph)
    protein_clusters = lapply(graph_decomposed, function(x) {
        vertices = sort(names(igraph::V(x)))
        list(ProteinNameOriginal = vertices,
             ProteinName = rep(paste(vertices, sep = ";", collapse = ";"),
                               times = length(vertices)))
    })
    protein_clusters = data.table::rbindlist(protein_clusters)
    protein_clusters
}

getIdenticalIndicator = function(pp_m) {
    res_m = matrix(0, nrow = ncol(pp_m), ncol = ncol(pp_m))
    colnames(res_m) = colnames(pp_m)
    for (i in seq_len(ncol(pp_m) - 1)) {
        for (j in (i + 1):ncol(pp_m)) {
            res_m[i, j] = all(pp_m[, i] == pp_m[, j])
        }
    }
    res_m = res_m + t(res_m)
    diag(res_m) = 1
    res_m
}
