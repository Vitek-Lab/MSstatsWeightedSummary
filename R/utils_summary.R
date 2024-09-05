#' Optimize a given criterion to estimate protein abundances
#'
#' @inheritParams getWeightedProteinSummary
#' @param weights data.table with weights.
#' @param use_shared if TRUE, shared peptides will be used.
#'
#' @export
#'
summarizeProteinsClusterSingleRun = function(feature_data, weights,
                                             norm, norm_parameter,
                                             use_shared) {
    run = unique(feature_data[, Run])
    feature_data = merge(feature_data[, list(ProteinName, PSM, Channel,
                                             log2IntensityNormalized)],
                         weights, by = c("ProteinName", "PSM"))
    proteins = unique(feature_data[["ProteinName"]])
    psms = unique(feature_data[["PSM"]])

    if (!use_shared) {
        feature_data[, IsUnique := data.table::uniqueN(ProteinName) == 1, by = "PSM"]
        feature_data = feature_data[(IsUnique)]
        feature_data[, Weight := 1]
        feature_data = unique(feature_data)
    }
    design = getProteinSummaryDesign(feature_data)
    design_matrix = design[["matrix"]]
    y = design[["response"]]
    is_non_missing = !is.na(y)
    optimization_problem = getProteinsOptimProblem(design_matrix[is_non_missing, ],
                                                   y[is_non_missing],
                                                   norm, norm_parameter,
                                                   proteins, psms)
    solution = CVXR::solve(optimization_problem)

    estimated_abundances = processProteinOptimSolution(solution,
                                                       optimization_problem,
                                                       feature_data,
                                                       design_matrix)
    estimated_abundances[, Run := run]
    estimated_abundances[, .(ProteinName, Run, Channel,
                             CenteredAbundance, Abundance)]
}

#' Get design matrix for protein-level summarization
#' @inheritParams getWeightedProteinSummary
#' @keywords internal
getProteinSummaryDesign = function(feature_data) {
    cols = c("PSM", "Channel", "log2IntensityNormalized")

    protein_intercepts = unique(feature_data[, .(ProteinName, PSM, Channel, Weight)])
    protein_intercepts = data.table::dcast(protein_intercepts,
                                           PSM + Channel ~ ProteinName,
                                           value.var = "Weight", fill = 0)

    feature_intercepts = unique(feature_data[, .(PSM, Channel, Present = 1)])
    feature_intercepts = data.table::dcast(feature_intercepts,
                                           PSM + Channel ~ PSM,
                                           value.var = "Present", fill = 0)

    channel_design = unique(feature_data[, .(ProteinName, PSM, Channel, Weight)])
    channel_design = data.table::dcast(channel_design,
                                       PSM + Channel ~ ProteinName + Channel,
                                       value.var = "Weight", fill = 0,
                                       sep = "__")

    intensities = unique(feature_data[, .(PSM, Channel, log2IntensityNormalized)])

    dm = merge(intensities, channel_design, by = c("PSM", "Channel"))
    dm = merge(dm, feature_intercepts, by = c("PSM", "Channel"))
    dm = merge(dm, protein_intercepts, by = c("PSM", "Channel"))
    dm[["intercept"]] = 1

    y = dm[["log2IntensityNormalized"]]
    dm = dm[, !(colnames(dm) %in% cols), with = FALSE]
    list(matrix = as.matrix(dm),
         response = y)
}

#' Get optimization problem for protein-level summaries
#' @param design_matrix design matrix for protein model optimization
#' @param y observed feature intensities
#' @inheritParams getWeightedProteinSummary
#' @param proteins vector of unique proteins
#' @param vector of unique psms
#' @keywords internal
getProteinsOptimProblem = function(design_matrix, y, norm,
                                   norm_parameter, proteins, psms) {
    num_params = ncol(design_matrix)

    protein_intercepts_condition = matrix(0, nrow = 1, ncol = num_params)
    protein_intercepts_condition[colnames(design_matrix) %in% proteins] = 1
    protein_profile_conditions = matrix(0, nrow = length(proteins),
                                        ncol = num_params)
    for (protein_id in seq_along(proteins)) {
        protein_cols = grepl(paste0(proteins[protein_id], "__"),
                             colnames(design_matrix))
        protein_profile_conditions[protein_id, protein_cols] = 1
    }

    feature_condition = matrix(0, nrow = 1, ncol = num_params)
    feature_condition[colnames(design_matrix) %in% psms] = 1

    constraints = rbind(protein_intercepts_condition,
                        protein_profile_conditions,
                        feature_condition)

    params = CVXR::Variable(num_params)
    if (norm == "p_norm") {
        objective = CVXR::Minimize(CVXR::p_norm(design_matrix %*% params - y,
                                                p = norm_parameter))
    } else if (norm == "Huber") {
        objective = CVXR::Minimize(sum(CVXR::huber(design_matrix %*% params - y,
                                                   M = norm_parameter)))
    }
    prob = CVXR::Problem(objective,
                         list(constraints %*% params == 0))
    prob
}

#' Get protein abundances from optimization result
#' @param solution solution generated by CVXR::solve
#' @param optimization_problem object returned by CVXR::Problem
#' @inheritParams getWeightedProteinSummary
#' @inheritParams getProteinsOptimProblem
#' @keywords internal
processProteinOptimSolution = function(solution, optimization_problem,
                                       feature_data, design_matrix) {
    num_channels = data.table::uniqueN(feature_data[, Channel])
    num_proteins = data.table::uniqueN(feature_data[, ProteinName])

    if (!(solution[["status"]] == "solver_error")) {
        param_names = colnames(design_matrix)
        estimated = solution$getValue(CVXR::variables(optimization_problem)[[1]])
        baseline = estimated[param_names == "intercept"]

        is_prot = param_names %in% unique(feature_data[, ProteinName])
        is_channel = grepl("__", param_names)

        protein_intercepts = data.table::data.table(
            ProteinName = param_names[is_prot],
            protein_intercept = estimated[is_prot]
        )

        protein_channel_names = stringr::str_split(param_names, "__")[is_channel]
        channels = estimated[is_channel]
        channels_dt = data.table::data.table(
            ProteinName = sapply(protein_channel_names,
                                 function(x) x[[1]]),
            Channel = sapply(protein_channel_names,
                             function(x) x[[2]]),
            ChannelValue = channels
        )

        result = merge(channels_dt, protein_intercepts, by = "ProteinName")
        result[, Intercept := baseline]
        result[, Abundance := ChannelValue + protein_intercept + Intercept]
        result[, CenteredAbundance := ChannelValue + protein_intercept]
        result[, list(ProteinName, Channel, CenteredAbundance, Abundance, Converged = TRUE)]
    } else {
        data.table::data.table(ProteinName = unique(feature_data[["ProteinName"]]),
                               Channel = unique(feature_data[["Channel"]]),
                               CenteredAbundance = NA_real_,
                               Abundance = NA_real_,
                               Converged = FALSE)
    }
}
