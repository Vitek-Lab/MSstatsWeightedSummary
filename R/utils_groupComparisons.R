#' Get groupComparison for different summarizations
#'
#' @param contrast_matrix contrast matrix
#' @param ... data.table object with summarization in MSstatsTMT format
#'
#' @return data.table
#'
#' @export
#'
getGroupComparisons = function(contrast_matrix, ...) {
    summaries_list = list(...)
    group_comparisons = lapply(summaries_list, function(x) {
        method = unique(x$Method)
        data.table::setnames(x, "ProteinName", "Protein", skip_absent = TRUE)
        group_comparison = MSstatsTMT::groupComparisonTMT(x, contrast_matrix)
        group_comparison$Method = method
        group_comparison
    })
    rbindlist(group_comparisons)
}