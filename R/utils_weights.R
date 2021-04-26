#' Get PSM-protein weights for summarization with shared peptides
#'
#' @param input data.table
#' @param method "all" for jointly estimating weights, "separate" for estimating
#' weights for each PSM separately
#' @param norm "p_norm" or "huber"
#' @param norm_parameter p for norm=="p_norm", M for norm=="Huber"
#'
#' @return data.table
#'
#' @export
#'
getAllWeights = function(input, method, norm, norm_parameter) {
    if (method == "all") {
        centered_summaries = unique(input[, .(ProteinName, Channel, Abundance)])
        centered_summaries[, AbundanceNew := Abundance - min(Abundance, na.rm = TRUE)]
        input = merge(input, centered_summaries[, .(ProteinName, Channel, AbundanceNew)],
                      by = c("ProteinName", "Channel"))
        input[, CenteredIntensity := log2IntensityNormalized - min(log2IntensityNormalized, na.rm = TRUE),
              by = c("ProteinName", "PSM")]
        input[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
        alphas_shared = .getWeightsCombined(input[!(IsUnique)], norm, norm_parameter)
        alphas_unique = unique(input[(IsUnique), .(ProteinName, PSM, Weight = 1)])
        alphas = rbind(alphas_shared, alphas_unique)
    } else {
        split_df = split(input, input$PSM)
        alphas = lapply(split_df, function(x) getPPWeights(x, TRUE, norm, norm_parameter))
        alphas = rbindlist(alphas)
    }
    alphas$Weight = ifelse(alphas$Weight < 0, 0, alphas$Weight)
    alphas
}


#' Get weights based on joint optimization (not per PSM)
#' @keywords internal
.getWeightsCombined = function(input, norm = "p_norm", norm_parameter = 1) {
    wide = data.table::dcast(input, CenteredIntensity + Channel ~ ProteinName + PSM,
                       value.var = "AbundanceNew", fill = 0)
    y_full = wide$CenteredIntensity
    x_full = as.matrix(wide[, -(1:2)])
    x_full = cbind(x_full, 1)

    cols = colnames(x_full)
    psms_cols = stringr::str_replace(cols, "BRD[0-9]_HUMAN_", "") # TODO: general proteins
    unique_psms = unique(psms_cols[-length(psms_cols)])
    n_params = ncol(x_full)
    n_conditions = length(unique_psms)

    constraint_matrix_full = matrix(0, nrow = n_conditions, ncol = n_params)
    for (i in seq_along(unique_psms)) {
        psm = unique_psms[i]
        constraint_matrix_full[i, psms_cols == psm] = 1
    }

    params_full = Variable(n_params)
    positive_matrix = diag(1, n_params - 1)
    positive_matrix = rbind(positive_matrix, 0)
    positive_matrix = cbind(positive_matrix, 0)
    constraints_full = list(positive_matrix %*% params_full >= rep(0, n_params),
                            constraint_matrix_full %*% params_full == rep(1, n_conditions))
    if (norm == "p_norm") {
        obj = CVXR::p_norm(x_full %*% params_full - y_full, norm_parameter)
    } else {
        obj = CVXR::huber(x_full %*% params_full - y_full, norm_parameter)
    }
    prob_con = CVXR::Problem(CVXR::Minimize(obj), constraints_full)
    sol_con = CVXR::solve(prob_con)
    alphas = as.vector(sol_con$getValue(variables(prob_con)[[1]]))

    result = data.table::data.table(
        ProteinName = stringr::str_extract(cols, "BRD[0-9]_HUMAN"), # TODO: general protein
        PSM = psms_cols,
        Weight = alphas
    )
    result[!is.na(ProteinName)]
}


#' Get weights for a single PSM
#' @keywords internal
getPPWeights = function(psm_data, intercept = TRUE, norm = "p_norm",
                        norm_parameter = 1) {
    wide = data.table::dcast(psm_data, log2IntensityNormalized + Channel ~ ProteinName,
                             value.var = "Abundance", fill = 0)
    num_proteins = ncol(wide) - 2

    if (num_proteins == 1) {
        return(data.table::data.table(PSM = unique(psm_data$PSM),
                                      ProteinName = unique(psm_data$ProteinName),
                                      Weight = 1,
                                      Optimization = "none"))
    } else {
        optimization_problem = .getWeightsOptimProblem(wide, intercept,
                                                       norm, norm_parameter)
        solution = CVXR::solve(optimization_problem)
        result = .processOptimSolution(solution, optimization_problem,
                                       num_proteins)
        result$PSM = unique(psm_data$PSM)
        result$ProteinName = colnames(wide)[3:ncol(wide)]
        result[, list(PSM, ProteinName, Weight, Optimization)]
    }
}


#' Get optimization problem for weights estimation
#' @keywords internal
.getWeightsOptimProblem = function(wide, intercept = TRUE, norm = "p_norm",
                                   norm_parameter = 1) {
    num_params = ncol(wide) - 2
    x_m = as.matrix(wide[, 3:ncol(wide), with = FALSE])
    x_m = apply(x_m, 2, function(x) x - mean(x))
    if (intercept) {
        params = CVXR::Variable(num_params + 1)
        x_m = cbind(x_m, 1)
    } else {
        params = CVXR::Variable(num_params)
    }
    y = wide$log2IntensityNormalized
    if (norm == "p_norm") {
        objective = CVXR::Minimize(CVXR::p_norm(x_m %*% params - y, p = norm_parameter))
    } else if (norm == "Huber") {
        objective = CVXR::Minimize(sum(CVXR::huber(x_m %*% params - y, M = norm_parameter)))
    }
    constraints = list(matrix(c(rep(1, num_params),
                                rep(0, as.integer(intercept))),
                              nrow = 1) %*% params == matrix(1, ncol = 1),
                       params >= rep(0, ncol(x_m)))
    prob = CVXR::Problem(objective, constraints)
    prob
}


#' Process output of optimization procedure
#' @keywords internal
.processOptimSolution = function(solution, problem, num_proteins) {
    if (solution$status == "optimal") {
        weights = solution$getValue(variables(problem)[[1]])
        weights = weights[1:num_proteins]
        status = "success"
    } else {
        weights = rep(1 / num_proteins, num_proteins)
        status = "failure"
    }
    data.table::data.table(Weight = weights,
                           Optimization = status)
}