#' Add missing proteins to PSM-level data based on a fasta database
#'
#' @param quantification_data MS data, preferably in `MSstats` or `MSstatsTMT` format.
#' @param fasta_db_path path to a fasta file that store protein sequences
#' @param protein_column name of a column with protein names.
#' @param peptide_column name of a column with peptide sequences.
#' @param n_cores number of cores that will be used while searching the database.
#' @param keep_unmodified if TRUE, a column that stores sequences of unmodified
#' peptides will be added to output data
#'
#' @importFrom data.table as.data.table data.table tstrsplit
#' @importFrom Biostrings vmatchPattern readAAStringSet
#' @importFrom IRanges elementNROWS
#'
#' @export
#'
adjustProteinAssignments = function(quantification_data, fasta_db_path,
                                    protein_column = "ProteinName",
                                    peptide_column = "PeptideSequence",
                                    n_cores = 1,
                                    keep_unmodified = FALSE) {
    quantification_data = data.table::as.data.table(quantification_data)
    protein_col = colnames(quantification_data) == protein_column
    quantification_data = quantification_data[, !protein_col, with = FALSE]

    all_peptide_sequences = unique(quantification_data[[peptide_column]])
    peptides_no_modifications = gsub("[^A-Z]", "", all_peptide_sequences)

    all_proteins = Biostrings::readAAStringSet(fasta_db_path)
    matching_proteins = data.table::rbindlist(
        parallel::mclapply(peptides_no_modifications, function(one_seq) {
            protein_indices = Biostrings::vmatchPattern(one_seq, all_proteins, fixed = TRUE)
            matching_proteins = names(protein_indices[IRanges::elementNROWS(protein_indices) > 0])
            data.table::data.table(PeptideSequenceUnmodified = one_seq,
                                   ProteinName = matching_proteins)
        }, mc.cores = n_cores)
    )

    matching_proteins[, ProteinName := data.table::tstrsplit(ProteinName, split = " ",
                                                             keep = 1)]
    matching_proteins = merge(data.table::data.table(
        PeptideSequence = all_peptide_sequences,
        PeptideSequenceUnmodified = peptides_no_modifications
    ),
    matching_proteins, by = "PeptideSequenceUnmodified", allow.cartesian = TRUE)
    matching_proteins = unique(matching_proteins)
    quantification_data = merge(quantification_data, matching_proteins,
                              by = "PeptideSequence", allow.cartesian = TRUE)
    if (!keep_unmodified) {
        unmodified_col = (colnames(quantification_data) == "PeptideSequenceUnmodified")
        quantification_data = quantification_data[, !unmodified_col, with = FALSE]
    }
    quantification_data[!is.na(ProteinName), ]
}
