#' Get a summary to start iterative optimization
#'
#' @param input data.table
#' @param method if "MSstats", MSstatsTMT summarization will be used. Otherwise,
#' summary will be based on "norm" and "norm_parameter" parameters.
#' @param use_discordant if TRUE, discordant peptides will be used
#' @param use_shared if TRUE, shared peptides will be used
#' @param norm "p_norm" or "huber"
#' @param norm_parameter p for norm=="p_norm", M for norm=="Huber"
#'
#' @return data.table
#' @import data.table
#'
#' @export
#'
getInitialSummary = function(input, method, use_discordant, use_shared, norm, norm_parameter) {
    if (method == "msstats") {
        result = getMSstatsSummary(input, use_shared, use_discordant)
    } else {
        weights_df = unique(input[, list(ProteinName, PSM, Weight)])
        result = getWeightedSummary(input, weights_df, use_discordant, norm,
                                    norm_parameter, use_shared, FALSE)
    }
    result = as.data.table(result)
    setnames(result, "Protein", "ProteinName", skip_absent = TRUE)
    result
}


#' Get MSstatsTMT-based protein-level summary
#' @keywords internal
getMSstatsSummary = function(input, use_shared, use_discordant) {
    input$Intensity = 2 ^ input$log2IntensityNormalized
    if (!use_discordant) {
        if (is.element("IsDiscordant", colnames(input))) {
            input = input[!(IsDiscordant)]
        }
    }
    input[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
    input_split = split(input, input$ProteinName)
    initial_summary = rbindlist(lapply(input_split, function(x) {
        x$Intensity = 2 ^ x$log2IntensityNormalized
        if (any(x$IsUnique)) {
            proteinSummarization(x[(IsUnique)], global_norm = FALSE,
                                 use_log_file = FALSE, verbose = FALSE)$ProteinLevelData
        } else {
            proteinSummarization(x, global_norm = FALSE,
                                 use_log_file = FALSE, verbose = FALSE)$ProteinLevelData
        }
    }))
    initial_summary
}


#' Optimize a give criterion to estimate protein abundances
#' @inheritParams getRobustSummary
#' @param weights_df data.table with weights.
#' @param use_shared if TRUE, shared peptides will be used.
#' @export
getWeightedSummary = function(input, weights_df, use_discordant,
                              norm, norm_parameter, use_shared, adaptive_huber) {
    input = merge(input[, list(ProteinName, PSM, Channel,
                               log2IntensityNormalized)],
                  weights_df, by = c("ProteinName", "PSM"))
    if (!use_discordant) {
        if (is.element("IsDiscordant", colnames(input))) {
            input = input[!(IsDiscordant)]
        }
    }
    if (!use_shared) {
        input[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
        input = input[(IsUnique)]
        input = unique(input)
    }
    design = .getDesign(input, "short")
    design_matrix = design[["matrix"]]
    y = design[["response"]]
    optimization_problem = .getProteinsOptimProblem(design_matrix, y,
                                                    norm, norm_parameter)
    solution = CVXR::solve(optimization_problem)

    if (norm == "Huber" & adaptive_huber) {
        M = norm_parameter
        prot_values = solution$getValue(CVXR::variables(optimization_problem)[[1]])
        raw_resid = design_matrix %*% prot_values - y
        huber_resid = ifelse(abs(raw_resid) < M, 0.5 * raw_resid ^ 2, M * (abs(raw_resid) - 0.5 * M))
        M = median(huber_resid)
    }

    estimated_abundances = .processProteinSolution(solution, optimization_problem,
                                                   input, design_matrix)
    if (norm == "Huber" & adaptive_huber) {
        attr(estimated_abundances, "M") = M
    }
    estimated_abundances
}

#' Get design matrix for protein-level summarization
#' @keywords internal
.getDesign = function(input, design_matrix) {
    cols = c("PSM", "Channel", "log2IntensityNormalized")
    dm = data.table::dcast(input, log2IntensityNormalized + PSM + Channel ~ ProteinName + Channel,
                           value.var = "Weight", fill = 0, sep = "__")
    dm$Present = 1
    dm = data.table::dcast(dm, ... ~ PSM, value.var = "Present", fill = 0)
    dm$intercept = 1
    y = dm$log2IntensityNormalized
    dm = as.matrix(dm[, !(colnames(dm) %in% cols), with = FALSE])
    list(matrix = dm,
         response = y)
}

#' Get optimization problem for protein-level summaries
#' @keywords internal
.getProteinsOptimProblem = function(design_matrix, y, norm = "p_norm",
                                    norm_parameter = 1) {
    num_params = ncol(design_matrix)
    params = CVXR::Variable(num_params)
    if (norm == "p_norm") {
        objective = CVXR::Minimize(CVXR::p_norm(design_matrix %*% params - y, p = norm_parameter))
    } else if (norm == "Huber") {
        objective = CVXR::Minimize(sum(CVXR::huber(design_matrix %*% params - y, M = norm_parameter)))
    }
    prob = CVXR::Problem(objective)
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

    param_names = colnames(design_matrix)
    result = data.table(label = param_names,
                        value = protein_values[, 1])
    protein_channel_names = stringr::str_split(result$label, "__")
    is_prot_channel = sapply(protein_channel_names, function(x) length(x) == 2)
    result = result[is_prot_channel]
    protein_channel_names = protein_channel_names[is_prot_channel]
    result[, ProteinName := sapply(protein_channel_names, function(x) x[1])]
    result[, Channel := sapply(protein_channel_names, function(x) x[2])]
    result[, Abundance := value + baseline]
    result[, list(ProteinName, Channel, Abundance)]
}
