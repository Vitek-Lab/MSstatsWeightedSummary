# MSstatsWeightedSummary: Weighted Protein-Level Summarization for Protein Clusters

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

TBA.
