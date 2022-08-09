#' Get PSM-protein weights for summarization with shared peptides
#'
#' @param input data.table
#' @param method "all" for jointly estimating weights, "separate" for estimating
#' weights for each PSM separately
#' @param norm "p_norm" or "huber"
#' @param norm_parameter p for norm=="p_norm", M for norm=="Huber"
#' @param equalize_protein_features if TRUE, all features in a protein group
#' that belong to a given protein will share weight parameter
#' @param weights_mode "contributions" for "sum to one" condition,
#' "probabilities" for only "non-negative" condition.
#'
#' @return data.table
#'
#' @export
#'
getAllWeights = function(input, method, norm, norm_parameter,
                         equalize_protein_features, weights_mode,
                         use_control, control_pattern) {
    if (!use_control) {
        input = input[!grepl(control_pattern, Condition)]
    }
    if (method == "all") {
        input[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
        alphas_shared = .getWeightsCombined(input, norm, norm_parameter,
                                            equalize_protein_features, weights_mode)
        # alphas_unique = unique(input[(IsUnique), .(ProteinName, PSM, Weight = 1)])
        # alphas = rbind(alphas_shared, alphas_unique)
        alphas  = alphas_shared
    } else {
        split_df = split(input, input$PSM)
        alphas = lapply(split_df, function(x) getPPWeights(x, TRUE, norm, norm_parameter,
                                                           weights_mode))
        alphas = rbindlist(alphas)
    }
    alphas$Weight = ifelse(alphas$Weight < 0, 0, alphas$Weight)
    alphas
}


#' Get weights based on joint optimization (not per PSM)
#' @keywords internal
.getWeightsCombined = function(input, norm = "p_norm", norm_parameter = 1,
                               equalize_protein_features = FALSE, weights_mode = "contributions") {
    intensities_tbl = unique(input[, .(PSM, Channel, log2IntensityNormalized)])
    intensities_tbl[, log2IntensityNormalized := log2IntensityNormalized - median(log2IntensityNormalized),
                    by = "PSM"]
    psms_intercept_tbl = unique(input[, .(PSM, Channel, Present = 1)])
    psms_intercept_tbl = dcast(psms_intercept_tbl,
                               PSM + Channel ~ PSM, value.var = "Present", fill = 0)
    colnames(psms_intercept_tbl)[3:ncol(psms_intercept_tbl)] = paste0("int_", 1:(ncol(psms_intercept_tbl) - 2))
    psms_protein_tbl = unique(input[, .(PSM, ProteinName, Channel, Abundance)])
    psms_protein_tbl[, Abundance := Abundance - median(unique(Abundance)),
                     by = "ProteinName"]
    psms_protein_tbl = dcast(psms_protein_tbl,
                             PSM + Channel ~ PSM + ProteinName,
                             value.var = "Abundance", fill = 0, sep = "___")
    wide = merge(
        merge(psms_intercept_tbl, intensities_tbl,
              by = c("PSM", "Channel"), all.x = T, all.y = T),
        psms_protein_tbl,
        by = c("PSM", "Channel"), all.x = T, all.y = T
    )

    y_full = wide$log2IntensityNormalized
    x_full = as.matrix(wide[, -(1:2), with = FALSE])
    x_full = x_full[, !(colnames(x_full) == "log2IntensityNormalized")]
    # x_full = cbind(intercept = 1, x_full)

    cols = colnames(x_full)
    cols_split = stringr::str_split(cols, "___")
    psms_cols = sapply(cols_split, function(x) x[1])
    protein_cols = sapply(cols_split, function(x) if (length(x) == 1) NA else x[-1])
    unique_psms = unique(psms_cols)
    unique_psms = unique_psms[!grepl("int", unique_psms)]
    n_params = ncol(x_full)
    n_conditions = length(unique_psms)

    constraint_matrix_full = matrix(0, nrow = n_conditions, ncol = n_params)
    for (i in seq_along(unique_psms)) {
        psm = unique_psms[i]
        constraint_matrix_full[i, psms_cols == psm] = 1
    }

    params_full = CVXR::Variable(n_params)
    n_intercepts = sum(grepl("int", cols))
    positive_matrix = diag(1, n_params - n_intercepts)
    positive_matrix = rbind(matrix(0, nrow = n_intercepts,
                                   ncol = ncol(positive_matrix)),
                            positive_matrix)
    positive_matrix = cbind(matrix(0, ncol = n_intercepts,
                                   nrow = nrow(positive_matrix)),
                            positive_matrix)

    protein_groups_dt = input[!(IsUnique),
                              .(AllProteins = paste(unique(ProteinName), sep = "__", collapse = "__")),
                              by = "PSM"]
    protein_groups = unique(protein_groups_dt$AllProteins)
    zeroes = matrix(0, nrow = 1, ncol = n_params)
    one_weight_constraint = do.call("rbind", unlist(unlist(lapply(protein_groups, function(protein_group) {
        proteins = stringr::str_split(protein_group, "__")[[1]]
        psms_group = protein_groups_dt[AllProteins == protein_group, unique(PSM)]
        lapply(proteins, function(protein) {
            first_psm = psms_group[[1]]
            lapply(psms_group[-1], function(psm) {
                x = zeroes
                x[, psms_cols %in% c(first_psm, psm) & protein_cols == protein] = c(-1, 1)
                x
            })
        })
    }), FALSE, FALSE), FALSE, FALSE))

    if (equalize_protein_features) {
        if (weights_mode == "contributions") {
            constraints_full = list(positive_matrix %*% params_full >= rep(0, n_params),
                                    constraint_matrix_full %*% params_full == c(rep(1, n_conditions),
                                                                                rep(0, nrow(constraint_matrix_full) - n_conditions)),
                                    one_weight_constraint %*% params_full == rep(0, nrow(one_weight_constraint)))
        } else {
            constraints_full = list(positive_matrix %*% params_full >= rep(0, n_params),
                                    positive_matrix %*% params_full <= rep(1, n_params),
                                    one_weight_constraint %*% params_full == rep(0, nrow(one_weight_constraint)))
        }
    } else {
        if (weights_mode == "contributions") {
            constraints_full = list(positive_matrix %*% params_full >= rep(0, n_params),
                                    constraint_matrix_full %*% params_full == c(rep(1, n_conditions),
                                                                                rep(0, nrow(constraint_matrix_full) - n_conditions)))
        } else {
            constraints_full = list(positive_matrix %*% params_full >= rep(0, n_params),
                                    positive_matrix %*% params_full <= rep(1, n_params))
        }

    }

    if (norm == "p_norm") {
        obj = CVXR::p_norm(x_full %*% params_full - y_full, norm_parameter)
    } else {
        obj = sum(CVXR::huber(x_full %*% params_full - y_full, norm_parameter))
    }
    # obj = CVXR::p_norm(x_full %*% params_full - y_full, 2)
    prob_con = CVXR::Problem(CVXR::Minimize(obj), constraints_full)
    sol_con = CVXR::solve(prob_con)
    alphas = as.vector(sol_con$getValue(CVXR::variables(prob_con)[[1]]))

    result = data.table::data.table(
        ProteinName = protein_cols, # TODO: general protein
        PSM = psms_cols,
        Weight = alphas
    )
    result[!is.na(ProteinName)]
}


#' Get weights for a single PSM
#' @keywords internal
getPPWeights = function(psm_data, intercept = FALSE, norm = "p_norm",
                        norm_parameter = 1, weights_mode = "contributions") {
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
                                                       norm, norm_parameter,
                                                       weights_mode)
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
.getWeightsOptimProblem = function(wide, intercept = FALSE, norm = "p_norm",
                                   norm_parameter = 1, weights_mode = "contributions") {
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
        objective = CVXR::Minimize(CVXR::p_norm(x_m %*% params - y, p = norm_parameter)) #+ CVXR::p_norm(params, 1))
    } else if (norm == "Huber") {
        objective = CVXR::Minimize(sum(CVXR::huber(x_m %*% params - y, M = norm_parameter)))# + CVXR::p_norm(params, 1))
    }
    if (weights_mode == "contributions") {
        constraints = list(matrix(c(rep(1, num_params),
                                    rep(0, as.integer(intercept))),
                                  nrow = 1) %*% params == matrix(1, ncol = 1),
                           params >= rep(0, ncol(x_m)))
    } else {
        constraints = list(params >= rep(0, ncol(x_m)),
                           params <= rep(1, ncol(x_m)))
    }
    prob = CVXR::Problem(objective, constraints)
    prob
}


#' Process output of optimization procedure
#' @keywords internal
.processOptimSolution = function(solution, problem, num_proteins) {
    if (solution$status == "optimal") {
        weights = solution$getValue(CVXR::variables(problem)[[1]])
        weights = weights[1:num_proteins]
        status = "success"
    } else {
        weights = rep(1 / num_proteins, num_proteins) # TODO: better solution
        status = "failure"
    }
    data.table::data.table(Weight = weights,
                           Optimization = status)
}

.getEqualProteinFeatureConstraint = function(input, cols, n_params) {
    pp = unique(input[, .(ProteinName, PSM)])
    pp[, Present := 1]
    pp[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
    pp = pp[!(IsUnique)]
    pp[, ProteinGroups := paste(unique(ProteinName), sep = ";", collapse = ";"),
       by = "PSM"]
    pp2 = unique(pp[, .(ProteinGroups, ProteinName, PSM, Present = 1)])
    pp3 = unique(pp2[order(ProteinGroups, PSM, ProteinName), .(ProteinGroups, PSM, ProteinName)])
    groups = split(pp3, list(pp3$ProteinGroups, pp3$ProteinName))
    groups = groups[sapply(groups, function(x) nrow(x) > 0)]

    constraints_by_group = vector("list", length(groups))
    for (i in seq_along(groups)) {
        group = groups[[i]]
        col_names_group = paste(group$ProteinName, group$PSM, sep = "_")
        reference = col_names_group[1]
        other_features = col_names_group[-1]
        group_result = matrix(0, nrow = length(other_features), ncol = n_params)
        for (j in seq_along(other_features)) {
            group_result[j, cols == reference] = 1
            group_result[j, cols == other_features[j]] = -1
        }
        constraints_by_group[[i]] = group_result
    }
    protein_constraint_matrix = do.call("rbind", constraints_by_group)
    colnames(protein_constraint_matrix) = cols
    protein_constraint_matrix
}

