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
#'
getAllWeights = function(input, norm, norm_parameter, weights_mode,
                         use_control, control_pattern) {
    input[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
    # input_stand = data.table::copy(input)
    # input_stand[, log2IntensityNormalize := log2IntensityNormalized - mean(log2IntensityNormalized),
    #             by = c("ProteinName", "PSM")]
    alphas_shared = .getWeightsCombined(input, norm, norm_parameter,
                                        weights_mode, use_control,
                                        control_pattern)
    # alphas_unique = unique(input[(IsUnique), .(ProteinName, PSM, Weight = 1)])
    # alphas = rbind(alphas_shared, alphas_unique)
    alphas  = alphas_shared
    alphas$Weight = ifelse(alphas$Weight < 0, 0, alphas$Weight)
    alphas
}


#' Get weights based on joint optimization (not per PSM)
#' @keywords internal
.getWeightsCombined = function(input, norm = "p_norm", norm_parameter = 1,
                               weights_mode = "contributions", use_control,
                               control_pattern) {
    # if (!use_control) {
    #     input = input[!grepl(control_pattern, Condition)]
    # }

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
