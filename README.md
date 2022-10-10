# Gene-Order-Workflow
A Nextflow workflow for gene order analysis of antimicrobial resistance (AMR) gene neighborhoods.

TO DO: Add introduction section explaining goals/utility of gene order based analysis.

## Workflow Summary 
* EXTRACTION: Identifies all AMR genes present in provided assemblies using Resistance Gene Identifier (RGI) annotations and extracts neighborhoods of a fixed window size N of genes upstream and downstream from a given focal AMR gene. 
* CLUSTERING: Derives siimilarity and distance matrices for each AMR gene, and applies three types of clustering algorithms (Hierarchical: UPGMA, Graph-based: MCL, Density-based: DBSCAN) to identify similarities and differences between AMR gene neighborhoods across genomes. 
* FILTERING (WIP): Collapses redundant (read: highly similar) gene neighborhoods prior to visualization for increased efficiency. Also filters low likelihood divergent genes prior to prediction phase. 
* PREDICTION (WIP): Creates gene-order based embeddings in the spirit of word embeddings and trains ML classifiers on the representations to predict candidate AMR genes from the set potential divergent AMR genes identified. Optional; must be enabled by the user when invoking workflow.  
* VISUALIZATION (WIP): Generates gene order visualizations, cluster visualizations, similarity/distance summary histograms, and interactive Dash Plotly dashboard for visual exploration of results. 

## Running the workflow 
The workflow requires user-provided paths to the user's RGI annotations, GBK files, results directory (for where final results will be stored), and BLAST files (if they exist; otherwise they will be automatically generated from the neighborhood FASTA files extracted). 

```bash
    nextflow run gene_neighborhoods_nf --rgi_path <RGI DIR> --gbk_path <GBK DIR> --output_path <RESULTS DIR> --blast_path <BLAST DIR> 
```

