#' Get PSM-protein weights for summarization with shared peptides
#'
#' @param input data.table
#' @param norm "p_norm" or "huber"
#' @param norm_parameter p for norm=="p_norm", M for norm=="Huber"
#' that belong to a given protein will share weight parameter
#' @param weights_mode "contributions" for "sum to one" condition,
#' "probabilities" for only "non-negative" condition.
#'
#' @return data.table
#'
#' @export
getAllWeights = function(input, norm = "p_norm", norm_parameter = 1,
                               weights_mode = "contributions", use_control = TRUE,
                               control_pattern = NULL) {
    input[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]

    if (!use_control) {
        input = input[!grepl(control_pattern, Condition)]
    }

    intensities_tbl = unique(input[, .(PSM, Channel, log2IntensityNormalized)])
    psms_intercept_tbl = unique(input[, .(PSM, Channel, Present = 1)])
    psms_intercept_tbl = data.table::dcast(psms_intercept_tbl,
                                           PSM + Channel ~ PSM, value.var = "Present", fill = 0)
    colnames(psms_intercept_tbl)[3:ncol(psms_intercept_tbl)] = paste0("int_", 1:(ncol(psms_intercept_tbl) - 2))
    psms_protein_tbl = unique(input[, .(PSM, ProteinName, Channel, Abundance)])
    psms_protein_tbl = dcast(psms_protein_tbl,
                             PSM + Channel ~ PSM + ProteinName,
                             value.var = "Abundance", fill = 0, sep = "___")
    wide = merge(
        merge(psms_intercept_tbl, intensities_tbl,
              by = c("PSM", "Channel"), all.x = T, all.y = T),
        psms_protein_tbl,
        by = c("PSM", "Channel"), all.x = T, all.y = T
    )
    wide = data.table::as.data.table(wide)

    y_full = wide$log2IntensityNormalized
    x_full = as.matrix(wide[, -(1:2), with = FALSE])
    x_full = x_full[, !(colnames(x_full) == "log2IntensityNormalized")]
    x_full = cbind(intercept = 1, x_full)

    cols = colnames(x_full)
    cols_split = stringr::str_split(cols, "___")
    psms_cols = sapply(cols_split, function(x) x[1])
    protein_cols = sapply(cols_split, function(x) if (length(x) == 1) NA else x[-1])
    unique_psms = unique(psms_cols)
    unique_psms = unique_psms[!grepl("int", unique_psms)]
    n_params = ncol(x_full)
    n_conditions = length(unique_psms)

    intercept_conditions = matrix(1, nrow = 1, ncol = n_params)
    intercept_conditions[!grepl("int_[0-9]+", colnames(x_full))] = 0

    constraint_matrix_full = matrix(0, nrow = n_conditions, ncol = n_params)
    for (i in seq_along(unique_psms)) {
        psm = unique_psms[i]
        constraint_matrix_full[i, psms_cols == psm] = 1
    }
    constraint_matrix_full = rbind(constraint_matrix_full,
                                   intercept_conditions)

    params_full = CVXR::Variable(n_params)
    n_intercepts = sum(grepl("int", cols))
    positive_matrix = diag(1, n_params - n_intercepts)
    positive_matrix = rbind(matrix(0, nrow = n_intercepts,
                                   ncol = ncol(positive_matrix)),
                            positive_matrix)
    positive_matrix = cbind(matrix(0, ncol = n_intercepts,
                                   nrow = nrow(positive_matrix)),
                            positive_matrix)

    if (weights_mode == "contributions") {
        constraints_full = list(positive_matrix %*% params_full >= rep(0, n_params),
                                constraint_matrix_full %*% params_full == c(rep(1, n_conditions),
                                                                            rep(0, nrow(constraint_matrix_full) - n_conditions)))
    } else {
        constraints_full = list(positive_matrix %*% params_full >= rep(0, n_params),
                                positive_matrix %*% params_full <= rep(1, n_params))
    }

    if (norm == "p_norm") {
        obj = CVXR::p_norm(x_full %*% params_full - y_full, norm_parameter)
    } else {
        obj = sum(CVXR::huber(x_full %*% params_full - y_full, norm_parameter))
    }
    prob_con = CVXR::Problem(CVXR::Minimize(obj), constraints_full)
    sol_con = CVXR::solve(prob_con)
    alphas = as.vector(sol_con$getValue(CVXR::variables(prob_con)[[1]]))

    result = data.table::data.table(
        ProteinName = protein_cols,
        PSM = psms_cols,
        Weight = alphas
    )
    result[!is.na(ProteinName)]
}
