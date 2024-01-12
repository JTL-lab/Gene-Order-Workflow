process PANAROO_RUN {
    //tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}"

    conda "bioconda::panaroo=1.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/panaroo:1.3.2--pyhdfd78af_0':
        'biocontainers/panaroo:1.3.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("results/*")                                      , emit: results
    tuple val(meta), path("results/core_gene_alignment.aln"), optional: true, emit: aln
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    panaroo \\
        $args \\
        -t $task.cpus \\
        -o results \\
        -i $gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panaroo: \$(echo \$(panaroo --version 2>&1) | sed 's/^.*panaroo //' ))
    END_VERSIONS
    """
}