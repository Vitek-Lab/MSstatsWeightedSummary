#' @importClassesFrom data.table data.table
setClassUnion("dtOrNULL", c("data.table", "NULL"))

#' Output of weighted summarization
#'
#' @slot FeatureLevelData feature-level (input) data.
#' @slot ProteinLevelData protein-level (summarized) output data.
#' @slot Weights a table of final peptide-protein Weights.
#' @slot FittedProfiles a table of observed and fitted PSM-level profiles.
#' @slot ConvergenceSummary a table with information about convergence for each Cluster and Run.
#' @slot FinalCriterionValues a table with values of the selected model-fitting criterion (p-norm or Huber loss) for converged models.
#' @slot WeightsHistory optional data.table of Weights from all iterations of fitting algorithm.
#' @slot ConvergenceHistory optional data.table with sums of absolute values of differences between Weights from consecutive iteration.
#' @slot ExperimentType "LF" or "TMT" depending on the type of input data.
#'
setClass("MSstatsWeightedSummary",
         slots = c(FeatureLevelData = "data.table",
                   ProteinLevelData = "data.table",
                   Weights = "data.table",
                   FittedProfiles = "data.table",
                   ConvergenceSummary = "data.table",
                   FinalCriterionValues = "data.table",
                   WeightsHistory = "dtOrNULL",
                   ConvergenceHistory = "dtOrNULL",
                   ExperimentType = "character"))

#' Extract feature-level data from MSstatsWeightedSummary object
#' @export
setGeneric("featureData",
           function(weighted_summary, proteins = NULL) standardGeneric("featureData"))
#' Extract feature-level data from MSstatsWeightedSummary object
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @param proteins optional character vector of proteins to extract. If NULL, all proteins
#' will be returned
#' @return data.table
setMethod("featureData", "MSstatsWeightedSummary",
          function(weighted_summary, proteins = NULL) {
              feature_level =  weighted_summary@FeatureLevelData
              if (!is.null(proteins)) {
                  feature_level = feature_level[ProteinName %in% proteins]
              }
              feature_level
          })


#' Extract protein-level data from MSstatsWeightedSummary object
#' @export
setGeneric("proteinData",
           function(weighted_summary, proteins = NULL) standardGeneric("proteinData"))
#' Extract protein-level data from MSstatsWeightedSummary object
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @param proteins optional character vector of proteins to extract. If NULL, all proteins
#' will be returned
#' @return data.table
setMethod("proteinData", "MSstatsWeightedSummary",
          function(weighted_summary, proteins = NULL) {
              protein_level =  weighted_summary@ProteinLevelData
              if (!is.null(proteins)) {
                  protein_level = protein_level[ProteinName %in% proteins]
              }
              protein_level
          })

#' Extract weights data from MSstatsWeightedSummary object
#' @export
setGeneric("featureWeights",
           function(weighted_summary, proteins = NULL, shared_only = TRUE)
               standardGeneric("featureWeights"))
#' Extract weights data from MSstatsWeightedSummary object
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @param proteins optional character vector of proteins to extract. If NULL, all proteins
#' will be returned
#' @param shared_only logical, output data.table will only include shared peptides
#' @return data.table
setMethod("featureWeights", "MSstatsWeightedSummary",
          function(weighted_summary, proteins = NULL, shared_only = TRUE) {
              weights = weighted_summary@Weights
              if (!is.null(proteins)) {
                  weights = weights[ProteinName %in% proteins]
              }
              if (shared_only) {
                  weights = weights[!(IsUnique)]
              }
              weights
          })

#' Extract convergence information from MSstatsWeightedSummary object
#' @export
setGeneric("convergenceSummary",
           function(weighted_summary) standardGeneric("convergenceSummary"))
#' Extract convergence information from MSstatsWeightedSummary object
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @return data.table
setMethod("convergenceSummary", "MSstatsWeightedSummary",
          function(weighted_summary) {
              weighted_summary@ConvergenceSummary
          })

#' Extract weights history from MSstatsWeightedSummary object
#' @export
setGeneric("weightsHistory",
           function(weighted_summary, shared_only = TRUE) standardGeneric("weightsHistory"))
#' Extract weights history from MSstatsWeightedSummary object
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @param shared_only logical, output data.table will only include shared peptides
#' @return data.table
setMethod("weightsHistory", "MSstatsWeightedSummary",
          function(weighted_summary, shared_only = TRUE) {
              weights_history = weighted_summary@WeightsHistory
              if (shared_only) {
                  weights_history = weights_history[!(IsUnique)]
              }
              weights_history
          })

#' Extract convergence history from MSstatsWeightedSummary object
#' @export
setGeneric("convergenceHistory",
           function(weighted_summary) standardGeneric("convergenceHistory"))
#' Extract convergence history from MSstatsWeightedSummary object
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @return data.table
setMethod("convergenceHistory", "MSstatsWeightedSummary",
          function(weighted_summary) {
              weighted_summary@ConvergenceHistory
          })

#' Extract cluster information from MSstatsWeightedSummary object
#' @export
setGeneric("proteinClusters",
           function(weighted_summary) standardGeneric("proteinClusters"))
#' Extract cluster information from MSstatsWeightedSummary object
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @return data.table
setMethod("proteinClusters", "MSstatsWeightedSummary",
          function(weighted_summary) {
              feature_data = weighted_summary@FeatureLevelData
              cluster_data = unique(feature_data[, .(Run, Cluster, ProteinName)])
              cluster_data
          })

#' Extract fitted PSM-level profiles from MSstatsWeightedSummary object
#' @export
setGeneric("fittedProfiles",
           function(weighted_summary) standardGeneric("fittedProfiles"))
#' Extract fitted PSM-level profiles from MSstatsWeightedSummary object
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @return data.table
setMethod("fittedProfiles", "MSstatsWeightedSummary",
          function(weighted_summary) {
              fitted_profiles = weighted_summary@FittedProfiles
              fitted_profiles
          })

#' Create input for MSstatsTMT::groupComparisonTMT function
#' @export
setGeneric("makeMSstatsTMTInput",
           function(weighted_summary, msstatstmt_output = NULL)
               standardGeneric("makeMSstatsTMTInput"))
#' Create input for MSstatsTMT::groupComparisonTMT function
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @param msstatstmt_output optional output of MSstatsTMT::proteinSummarization function
#' @return list
setMethod("makeMSstatsTMTInput", "MSstatsWeightedSummary",
          function(weighted_summary, msstatstmt_output = NULL) {
              feature_data = weighted_summary@FeatureLevelData
              protein_data = weighted_summary@ProteinLevelData

              if (is.null(msstatstmt_output)) {
                  list(FeatureLevelData = feature_data,
                       ProteinLevelData = protein_data)
              } else {
                  mstmt_feature = msstatstmt_output[["FeatureLevelData"]]
                  mstmt_protein = msstatstmt_output[["ProteinLevelData"]]
                  list(FeatureLevelData = rbind(feature_data,
                                                mstmt_feature,
                                                use.names = TRUE, fill = TRUE),
                       ProteinLevelData = rbind(protein_data,
                                                mstmt_protein,
                                                use.names = TRUE, fill = TRUE))
              }
          })

#' Create input for MSstats::groupComparison function
#' @export
setGeneric("makeMSstatsInput",
           function(weighted_summary, msstats_output = NULL)
               standardGeneric("makeMSstatsInput"))
#' Create input for MSstats::groupComparison function
#' @param weighted_summary output of the getWeightedProteinSummary function
#' @param msstatstmt_output optional output of MSstats::dataProcess function
#' @return list
setMethod("makeMSstatsInput", "MSstatsWeightedSummary",
          function(weighted_summary, msstats_output = NULL) {
              feature_data = data.table::copy(weighted_summary@FeatureLevelData)
              protein_data = data.table::copy(weighted_summary@ProteinLevelData)

              if (!is.null(msstats_output)) {
                  mst_feature = msstats_output[["FeatureLevelData"]]
                  mst_protein = msstats_output[["ProteinLevelData"]]
                  feature_data_conv = rbind(feature_data,
                                            mst_feature,
                                            use.names = TRUE, fill = TRUE)
                  protein_data = rbind(protein_data,
                                       mst_protein,
                                       use.names = TRUE, fill = TRUE)
              }
              if (weighted_summary@ExperimentType == "LF") {
                  feature_data_conv = feature_data[, list(
                      PROTEIN = ProteinName,
                      PEPTIDE = paste(PeptideSequence, Charge, sep = "_"),
                      TRANSITION = paste(FragmentIon, ProductCharge, sep = "_"),
                      FEATURE = PSM,
                      LABEL = IsotopeLabelType,
                      GROUP = Condition,
                      RUN = Run,
                      SUBJECT = BioReplicate,
                      FRACTION = Fraction,
                      originalRUN = Run,
                      censored = FALSE,
                      INTENSITY = Intensity,
                      ABUNDANCE = log2IntensityNormalized,
                      newABUNDANCE = log2IntensityNormalized,
                      predicted = NA_real_,
                      remove = FALSE)]
                  cols = c("RUN", "Protein", "LogIntensities", "originalRUN",
                           "GROUP", "SUBJECT", "more50missing", "NumMeasuredFeature")
                  setnames(protein_data,
                           c("Run", "Abundance", "Condition", "BioReplicate"),
                           c("originalRUN", "LogIntensities", "GROUP", "SUBJECT"))
                  protein_data[, RUN := originalRUN]
                  protein_data[, more50missing := FALSE]
                  num_features = feature_data[!is.na(log2IntensityNormalized),
                                              list(NumMeasuredFeature = data.table::uniqueN(PSM)),
                                              by = c("ProteinName", "Run")]
                  group_meas = feature_data[!is.na(log2IntensityNormalized),
                                            list(TotalGroupMeasurements = data.table::uniqueN(PSM)),
                                            by = c("ProteinName", "Condition")]
                  protein_data = merge(protein_data, num_features,
                                       by.x = c("Protein", "originalRUN"),
                                       by.y = c("ProteinName", "Run"),
                                       all.x = TRUE, all.y = TRUE, sort = FALSE)
                  protein_data = merge(protein_data, group_meas,
                                       by.x = c("Protein", "GROUP"),
                                       by.y = c("ProteinName", "Condition"),
                                       sort = FALSE)
                  protein_data[, MissingPercentage := 0.0]
                  protein_data[, NumImputedFeature := 0L]
              } else {
                  feature_data_conv = feature_data
              }
              list(FeatureLevelData = feature_data_conv,
                   ProteinLevelData = protein_data)
          })
