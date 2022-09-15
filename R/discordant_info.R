#' Cluster peptide profiles
#'
#' @param input data.table in MSstatsTMT format
#' @param use_control if TRUE, control group channels will be used in clustering
#' @param control_pattern regular expression that identifies control conditions
#'
#' @return data.table
#'
#' @export
#'
getDiscordantInfo = function(input, use_control = TRUE, control_pattern = "DMSO") {
    input_by_prot = split(input, input$ProteinName)
    results = lapply(input_by_prot, function(x) {
        profiles = unique(x[, .(PSM, Channel, Condition, log2IntensityNormalized)])
        if (!use_control) {
            profiles = profiles[!grepl(control_pattern, Condition)]
        }
        profiles = profiles[, .(PSM, Channel, log2IntensityNormalized)]
        wide = dcast(profiles, PSM ~ Channel, value.var = "log2IntensityNormalized",
                     fill = 0)
        psms = wide$PSM
        wide = wide[, -1]
        wide = as.matrix(wide)
        dists_cor = proxy::dist(wide, method = function(x, y) 1 - cor(x, y))
        clusters = try(NbClust::NbClust(data = wide, diss = dists_cor,
                                        min.nc = 1,
                                        max.nc = min(8, uniqueN(psms)),
                                        method = "average", distance = NULL,
                                        index = "gap"),
                       silent = TRUE)
        if (inherits(clusters, "try-error")) {
            result = data.table(PSM = psms,
                                PatternCluster = 1,
                                ClusterStatus = "failure")
        } else {
            result = data.table(PSM = psms,
                                PatternCluster = clusters$Best.partition,
                                ClusterStatus = "success")
        }
        merge(x, result, by = "PSM")
    })
    results = rbindlist(results)
    results$IsDiscordant = results$PatternCluster != 1
    results
}
