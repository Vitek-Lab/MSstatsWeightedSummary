#' @keywords internal
imputeTMT = function(feature_data) {
    censored_cutoff = 0.999

    if (any(is.na(feature_data[["log2IntensityNormalized"]]))) {
        quantiles = feature_data[!is.na(log2IntensityNormalized) & Intensity > 1,
                                 .(prob = c(0.01, 0.25, 0.5, 0.75,
                                            censored_cutoff),
                                   quant = quantile(log2IntensityNormalized,
                                                    prob = c(0.01, 0.25, 0.5, 0.75,
                                                             censored_cutoff),
                                                    na.rm = TRUE)), by = c("Run", "Cluster")]
        quantile_stats = quantiles[, .(iqr = quant[prob == 0.75] - quant[prob == 0.25],
                                       q25 = quant[prob == 0.25],
                                       mult = quant[prob == censored_cutoff] - quant[prob == 0.75]),
                                   by = c("Run", "Cluster")]
        quantile_stats[, cutoff_lower := q25 - mult]

        feature_data = merge(feature_data, quantile_stats[, .(Run, Cluster, cutoff_lower)], by = c("Run", "Cluster"))
        feature_data[, censored := (!is.na(log2IntensityNormalized) & log2IntensityNormalized < cutoff_lower) | is.na(log2IntensityNormalized)]

        feature_data[, nonmissing := !is.na(log2IntensityNormalized)]
        feature_data[, n_obs := sum(nonmissing), by = c("Run", "PSM")]
        feature_data[, nonmissing := ifelse(n_obs <= 1, FALSE, nonmissing)]
        feature_data[, n_obs_run := sum(nonmissing), by = c("Run", "Channel")]
        feature_data[, nonmissing_all := !is.na(log2IntensityNormalized)]
        feature_data[, total_features := uniqueN(PSM), by = c("Run")]
        feature_data[, n_obs := .N, by = c("Run", "PSM")]
        feature_data[, nonmissing_all := ifelse(total_features > 1 & n_obs <= 1,
                                                FALSE, nonmissing_all)]
        grouping_vars = c("Run", "PSM") # Per protein
        # grouping_vars = c("Run", "PSM") # Per cluster
        feature_data[n_obs > 1 & n_obs_run > 0,
                     ABUNDANCE_cut := getMin(log2IntensityNormalized, nonmissing_all),
                     by = grouping_vars]
        feature_data[, any_censored := any(censored & n_obs > 1 & n_obs_run > 0),
                     by = c("Run")]
        feature_data[, log2IntensityOrig := log2IntensityNormalized]
        feature_data[, log2IntensityNormalized := ifelse(!nonmissing_all & censored & is.finite(ABUNDANCE_cut) & any_censored,
                                                         ABUNDANCE_cut, log2IntensityNormalized)]
        feature_data[, cen := ifelse(censored, 0, 1)]
    }
    data_by_run_cluster = split(feature_data, feature_data[, .(Run, Cluster)])
    data_by_run_cluster = lapply(data_by_run_cluster, imputeSingleRunCluster)
    data_by_run_cluster = data.table::rbindlist(data_by_run_cluster, fill = TRUE)
    data_by_run_cluster[, log2IntensityNormalized := ifelse(!is.na(log2IntensityNormalized) & log2IntensityNormalized < 0,
                                                            NA_real_, log2IntensityNormalized)]
    data_by_run_cluster
}

#' @keywords internal
#' @importFrom survival Surv
imputeSingleRunCluster = function(feature_data) {
    if (any(is.na(feature_data$log2IntensityNormalized))) {
        feature_data[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
        feature_wide = dcast(unique(feature_data[, .(PSM, Channel, ProteinName, IsUnique,
                                                     log2IntensityNormalized, censored, cen, Present = 1)]),
                             PSM + Channel + log2IntensityNormalized + IsUnique + censored + cen ~ ProteinName,
                             value.var = "Present", fill = 0)
        feature_wide[, PSM := factor(PSM)]
        feature_wide[, Channel := factor(Channel)]

        if (uniqueN(feature_wide$PSM) > 1) {
            formula_chr = "Surv(log2IntensityNormalized, cen, type = 'left') ~ PSM + Channel"
        } else {
            formula_chr = "Surv(log2IntensityNormalized, cen, type = 'left') ~ Channel"
        }

        fit = survival::survreg(as.formula(formula_chr), data = feature_wide, dist = "gaussian")
        if (fit$df.residual >= 0) {
            feature_wide[, Imputed := predict(fit, newdata = feature_wide)]
        } else {
            feature_wide[, Imputed := NA_real_]
        }

        feature_data = merge(feature_data, unique(feature_wide[, .(PSM, Channel, Imputed)]),
                             by = c("PSM", "Channel"), all.x = T)
        feature_data[, log2IntensityNormalized := ifelse(is.na(log2IntensityNormalized), Imputed, log2IntensityNormalized)]
    }
    feature_data
}

#' @keywords internal
getMin = function(abundance, nonmissing) {
    0.99 * min(abundance[nonmissing], na.rm = TRUE)
}
