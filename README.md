# MSstatsWeightedSummary: Weighted Protein-Level Summarization for Protein Clusters

Protein-level summarization of feature (PSM)-level data is an important step in proteomics data analysis methods, including the MSstats workflow. MSstats and MSstatsTMT packages summarize features based on unique peptides. This package implements a statistical method to jointly summarize clusters of proteins that share features using weights specific to peptide-protein matches.

#### Installation:


```
if (!require(devtools)) {
    install.packages("devtools")
}
devtools::install_github("Vitek-Lab/MSstatsWeightedSummary")
```


#### Getting help:

[Vignette online](https://vitek-lab.github.io/MSstatsWeightedSummary/articles/shared_workflow.html)

```
?getWeightedProteinSummary
?normalizeSharedPeptides
?addClusterMembership
vignette(package = "MSstatsWeightedSummary")
```

#### Major functionalities:

  - `normalizeSharedPeptides`: feature-level data MSstatsTMT normalization that accounts for shared peptides
  - `addClusterMembership`: identification of clusters of proteins that share peptides 
    (connected components of peptide-protein graph)
  - `getWeightedSummary`: protein-level summarization of feature-level data. Performed separately for each
    cluster of proteins and each run. Input and output compatible with MSstatsTMT workflow for analysis
    of mass spectrometry proteomics data with isobaric labeling.
    
    
#### How to cite

This package implements a statistical model introduced in the article [Relative quantification of proteins and post-translational modifications in proteomic experiments with shared peptides: a weight-based approach](https://doi.org/10.1093/bioinformatics/btaf046). 
The bibtex entry can be found by calling:

```
citation(package = "MSstatsWeightedSummary")
```
