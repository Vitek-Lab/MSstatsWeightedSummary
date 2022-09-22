#' Optimize a given criterion to estimate protein abundances
#' @inheritParams getRobustSummary
#' @param weights_df data.table with weights.
#' @param use_shared if TRUE, shared peptides will be used.
#' @export
getWeightedSummary = function(input, weights_df,
                              norm, norm_parameter,
                              use_shared) {
    input = merge(input[, list(ProteinName, PSM, Channel,
                               log2IntensityNormalized)],
                  weights_df, by = c("ProteinName", "PSM"))
    proteins = unique(input[["ProteinName"]])
    psms = unique(input[["PSM"]])

    if (!use_shared) {
        input[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
        input[, HasUnique := any(IsUnique), by = "ProteinName"]
        input = input[(IsUnique) | !(HasUnique)]
        input = unique(input)
    }
    design = .getDesign(input)
    design_matrix = design[["matrix"]]
    y = design[["response"]]
    optimization_problem = .getProteinsOptimProblem(design_matrix, y,
                                                    norm, norm_parameter,
                                                    proteins, psms)
    solution = CVXR::solve(optimization_problem)

    estimated_abundances = .processProteinSolution(solution, optimization_problem,
                                                   input, design_matrix)
    estimated_abundances
}

#' Get design matrix for protein-level summarization
#' @keywords internal
.getDesign = function(input) {
    cols = c("PSM", "Channel", "log2IntensityNormalized")

    protein_intercepts = unique(input[, .(ProteinName, PSM, Channel, Weight)])
    protein_intercepts = data.table::dcast(protein_intercepts,
                                           PSM + Channel ~ ProteinName,
                                           value.var = "Weight", fill = 0)

    feature_intercepts = unique(input[, .(PSM, Channel, Present = 1)])
    feature_intercepts = data.table::dcast(feature_intercepts,
                                           PSM + Channel ~ PSM,
                                           value.var = "Present", fill = 0)

    channel_design = unique(input[, .(ProteinName, PSM, Channel, Weight)])
    channel_design = data.table::dcast(channel_design,
                                       PSM + Channel ~ ProteinName + Channel,
                                       value.var = "Weight", fill = 0,
                                       sep = "__")

    intensities = unique(input[, .(PSM, Channel, log2IntensityNormalized)])

    dm = merge(intensities, channel_design, by = c("PSM", "Channel"))
    dm = merge(dm, feature_intercepts, by = c("PSM", "Channel"))
    dm = merge(dm, protein_intercepts, by = c("PSM", "Channel"))
    dm[["intercept"]] = 1

    proteins = unique(input[["ProteinName"]])
    y = dm[["log2IntensityNormalized"]]
    dm = dm[, !(colnames(dm) %in% cols), with = FALSE]
    list(matrix = as.matrix(dm),
         response = y)
}

#' Get optimization problem for protein-level summaries
#' @keywords internal
.getProteinsOptimProblem = function(design_matrix, y, norm,
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

    conditions = rbind(protein_intercepts_condition,
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
                         list(conditions %*% params == 0))
    prob
}

#' Get protein abundances from optimization result
#' @keywords internal
.processProteinSolution = function(solution, optimization_problem,
                                   input, design_matrix) {
    num_channels = uniqueN(input$Channel)
    num_proteins = uniqueN(input$ProteinName)
    protein_values = solution$getValue(CVXR::variables(optimization_problem)[[1]])
    baseline = protein_values[length(protein_values)]
    protein_intercepts = protein_values[colnames(design_matrix) %in% unique(input$ProteinName)]

    param_names = colnames(design_matrix)
    result = data.table(label = param_names,
                        value = protein_values[, 1])
    protein_channel_names = stringr::str_split(result$label, "__")
    is_prot_channel = sapply(protein_channel_names, function(x) length(x) == 2)

    protein_intercepts = result[label %in% unique(input[["ProteinName"]])]
    setnames(protein_intercepts, c("ProteinName", "intercept"))

    result = result[is_prot_channel]
    protein_channel_names = protein_channel_names[is_prot_channel]
    result[, ProteinName := sapply(protein_channel_names, function(x) x[1])]
    result[, Channel := sapply(protein_channel_names, function(x) x[2])]
    result = merge(result, protein_intercepts, by = "ProteinName")
    result[, global_intercept := baseline]
    result[, Abundance := value + intercept]
    result[, list(ProteinName, Channel, Abundance, global_intercept)]
}
