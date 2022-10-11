# Gene-Order-Workflow
A Nextflow workflow for gene order analysis of antimicrobial resistance (AMR) gene neighborhoods.

TO DO: Add introduction section briefly explaining goals/utility of gene order based analysis.

## Workflow Summary 
* EXTRACTION: Identifies all AMR genes present in provided assemblies using Resistance Gene Identifier (RGI) annotations and extracts neighborhoods of a fixed window size N of genes upstream and downstream from a given focal AMR gene. 
* CLUSTERING: Derives siimilarity and distance matrices for each AMR gene, and applies three types of clustering algorithms (Hierarchical: UPGMA, Graph-based: MCL, Density-based: DBSCAN) to identify similarities and differences between AMR gene neighborhoods across genomes. 
* FILTERING (WIP): Collapses redundant (read: highly similar) gene neighborhoods prior to visualization for increased efficiency. Also filters low likelihood divergent genes prior to prediction phase. 
* PREDICTION (WIP): Creates gene-order based embeddings in the spirit of word embeddings and trains ML classifiers on the representations to predict candidate AMR genes from the set potential divergent AMR genes identified. Optional; must be enabled by the user when invoking workflow.  
* VISUALIZATION (WIP): Generates gene order visualizations, cluster visualizations, similarity/distance summary histograms, and interactive Dash Plotly dashboard for visual exploration of results. 

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Start running your own analysis!

   ```console
   nextflow run main.nf -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --rgi_path <RGI DIR> --gbk_path <GBK DIR> --output_path <OUTPUT DIR>
   ```

## Credits

nf-core/geneorderanalysis was originally written by Julia Lewandowski.


## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
