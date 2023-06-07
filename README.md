# Gene-Order-Workflow
A Nextflow workflow for bacterial gene order analysis, with outputs easily explorable through its partner visualization application [Coeus](https://github.com/JTL-lab/Coeus).

## Workflow Summary 
* EXTRACTION: Identifies all genes of interest present in provided assemblies and extracts neighborhoods consisting of `num_neighbors` genes upstream and downstream from a given focal gene. 
* CLUSTERING: Derives similarity and distance matrices for each AMR gene, and applies three types of clustering algorithms (Hierarchical: UPGMA, Graph-based: MCL, Density-based: DBSCAN) to identify similarities and differences between gene neighborhoods across genomes. 
* FILTERING: Collapses identical gene neighborhoods prior to visualization for increased efficiency. It retains one surrogate neighborhood instead of multiple identical ones, and outputs a textfile for each gene that had neighborhoods collapsed to indicate which genome is standing in for which other genomes.
* VISUALIZATION: Pre-computes gene order visualizations, cluster visualizations, and similarity/distance summary histograms. Clustering visualizations can also be dynamically generated with [Coeus](https://github.com/JTL-lab/Coeus) using similarity and distance matrices calculated through this workflow.
* (TBD) PREDICTION: If analyzing antimicrobial resistance (AMR) genes using annotations from AMR detection software, create gene-order based embeddings in the spirit of word embeddings and train ML classifiers on the representations to predict candidate AMR genes. Optional; must be enabled by the user when invoking workflow.  


## Data
Gene-Order-Workflow can be run on either: 
* Genbank files and their corresponding assemblies.
* Genbank files, their corresponding assemblies, and annotations from external software (e.g. the [Resistance Gene Identifier](https://github.com/arpcard/rgi)). 

#### Required Files
* INPUT FILE (.txt)
    * If using Genbank files only: input file should be a textfile with one gene name per row for every gene you want to analyze.
    ```
    geneA 
    geneB
    ...
    geneX
    ```
    * If using Genbank files and annotations: input file should be a textfile indicating the column names (values between '<>' should be replaced with your column names) in your annotation file that correspond to the `Contig` identifier, `Gene_Start`, `Gene_Stop`, and `Gene_Name` for your data. An example that can be used with [RGI](https://github.com/arpcard/rgi) can be found in `sample_data/rgi_input.txt`.
    ```
    <contig_col> = Contig
    <gene_name_col> = Gene_Name
    <start_col> = Gene_Start
    <stop_col> = Gene_End
    ```
* ASSEMBLIES ( .faa | .fa | .fna ) in a single directory. Should be annotated and share the same locus tags as those found in the Genbank files. 
* GENBANK FILES ( .gbk | .gb ) in a single directory. Should correspond to the assembly files. 

#### Optional Files
* ANNOTATIONS (e.g. RGI textfiles) in a single directory. 

## Workflow Parameters

| Parameter       | Type     | Description                                                                                                                                                                                                                                   |
|-----------------|----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| input_file_path | Required | Path to input textfile containing either a) gene names to extract if providing Genbank files and assemblies or b) column names of required columns if additionally providing annotations (see examples above).                                |
| assembly_path   | Required | Path to directory containing assembly files (.faa, .fa, .fna).                                                                                                                                                                                |
| gbk_path        | Required | Path to directory containing Genbank files (.gbk, .gb).                                                                                                                                                                                       |
| extract_path    | Optional | Path to annotation textfiles (.txt).                                                                                                                                                                                                          |
| num_neighbors   | Optional | Neighborhood size to extract. Should be an even number N, such that N/2 neighbors upstream and N/2 neighbors downstream will be analyzed. Default: 10.                                                                                        |
| percent_cutoff  | Optional | Cutoff percentage of genomes a gene should be present within to be included in extraction and subsequent analysis. Should a float between 0 and 1 (e.g., 0.25 means only genes present in a minimum of 25% of genomes are kept). Default: 0.25. |
| inflation       | Optional | Inflation hyperparameter value for Markov Clustering Algorithm. See the [algorithm documentation](https://markov-clustering.readthedocs.io/en/latest/readme.html) for details. Default: 2.                                                    |
| epsilon         | Optional | Epsilon hyperparameter value for DBSCAN clustering. See the [algorithm documentation](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html) for details. Default: 0.5.                                               | 
| minpts | Optional | Minpts hyperparameter value for DBSCAN clustering. See the [algorithm documentation](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html) for details. Default: 5.                                                  | 
| outdir | Optional | Path to output directory. Default: 'results' within repository. |                                                                                                                                                                                | 

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility (you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles)).

3. Start running your own analysis!
    ```console
    nextflow run main.nf -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> \
    --input_file_path <path_to_input_file> \
    --assembly_path <path_to_assembly_dir> \
    --extract_path <path_to_annotation_dir> \ 
    --gbk_path <path_to_genbank_dir> \
    --num_neighbors <int_val> --percent_cutoff <float_val>
    ```
   
Please ensure you've formatted your input file correctly for the use case you need (see Data section above).

Some general notes regarding running on HPC environments: 
* Please check nf-core/configs to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use -profile <institute> in your command. This will enable either docker or singularity and set the appropriate execution settings for your local compute environment.
* If you are using Singularity, then the pipeline will auto-detect this and attempt to download the Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Singularity images directly due to timeout or network issues then please use the --singularity_pull_docker_container parameter to pull and convert the Docker image instead. 
* When running Nextflow on HPC environments with the Slurm executor, some have reported persistent SIGBUS errors. If describes you, you may find it helpful to consult this [suggested fix](https://github.com/nextflow-io/nextflow/issues/842#issuecomment-567119760).

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

