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
    if (method == "mstats") {
        result = getMSstatsSummary(input, use_shared, use_discordant)
    } else {
        result = getOptimSummary(input, "short", use_discordant, norm,
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
        if (is.element(colnames(input), "IsDiscordant")) {
            input = input[!(IsDiscordant)]
        }
    }
    if (use_shared) {
        result = lapply(split(input, input[, .(ProteinName)]),
                        function(x) MSstatsTMT::proteinSummarization(
                            unique(x[, .(ProteinName, PeptideSequence, Charge, PSM, Channel,
                                         BioReplicate, Condition, Run, Mixture, TechRepMixture, Intensity)]),
                            global_norm = FALSE))
        result = rbindlist(result)
    } else {
        input[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
        input = input[(IsUnique)]
        input = unique(input[, .(ProteinName, PeptideSequence, Charge, PSM, Channel,
                                 BioReplicate, Condition, Run, Mixture, TechRepMixture, Intensity)])
        result = MSstatsTMT::proteinSummarization(input, global_norm = FALSE)
    }
    result
}


#' Optimize a give criterion to estimate protein abundances
#' @keywords internal
getOptimSummary = function(input, design_type, use_discordant,
                           norm, norm_parameter, use_shared, adaptive_huber) {
    if (!use_discordant) {
        if (is.element(colnames(input), "IsDiscordant")) {
            input = input[!(IsDiscordant)]
        }
    }
    if (!use_shared) {
        input[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
        input = input[(IsUnique)]
        input = unique(input)
    }
    design = .getDesign(input, design_type)
    design_matrix = design[["matrix"]]
    y = design[["response"]]
    optimization_problem = .getProteinsOptimProblem(design_matrix, y,
                                                    norm, norm_parameter)
    solution = CVXR::solve(optimization_problem)

    if (norm == "Huber" & adaptive_huber) {
        M = norm_parameter
        prot_values = solution$getValue(variables(optimization_problem)[[1]])
        raw_resid = design_matrix %*% prot_values - y
        huber_resid = ifelse(abs(raw_resid) < M, 0.5 * raw_resid ^ 2, M * (abs(raw_resid) - 0.5 * M))
        M = median(huber_resid)
    }

    estimated_abundances = .processProteinSolution(solution, optimization_problem,
                                                   input)
    if (norm == "Huber" & adaptive_huber) {
        attr(estimated_abundances, "M") = M
    }
    estimated_abundances
}

#' Get design matrix for protein-level summarization
#' @keywords internal
.getDesign = function(input, design_matrix) {
    cols = c("PSM", "Channel", "log2IntensityNormalized")
    if (design_matrix == "short") {
        dm = data.table::dcast(input, log2IntensityNormalized + PSM + Channel ~ ProteinName + Channel,
                   value.var = "Weight", fill = 0)
        dm$Present = 1
        dm = data.table::dcast(dm, ... ~ PSM, value.var = "Present", fill = 0)
        dm$intercept = 1
        y = dm$log2IntensityNormalized
        dm = as.matrix(dm[, !(colnames(dm) %in% cols), with = FALSE])
    } else {
        stop("Not implemented yet")
    }
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
                                   input) {
    num_channels = uniqueN(input$Channel)
    num_proteins = uniqueN(input$ProteinName)
    protein_values = solution$getValue(variables(optimization_problem)[[1]])
    baseline = protein_values[length(protein_values)]
    protein_values = baseline + protein_values[1:(num_channels * num_proteins)]
    data.table(ProteinName = rep(unique(input$ProteinName), each = num_channels),
               Channel = rep(unique(input$Channel), times = num_proteins),
               Abundance = protein_values)
}
