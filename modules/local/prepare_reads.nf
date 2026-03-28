process PREPARE_READS {
    tag "$meta.id"
    label 'process_low'

    // Uses only standard shell tools (gzip/zcat) — available in any container
    // Use the python container already pulled for other steps
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.fastq.gz"), emit: fastq_gz
    tuple val(meta), path("${meta.id}.fastq"),    emit: fastq
    path "versions.yml",                           emit: versions

    script:
    def is_gz = fastq.name.endsWith('.gz')
    """
    if ${is_gz}; then
        ln -s ${fastq} ${meta.id}.fastq.gz
        zcat ${fastq} > ${meta.id}.fastq
    else
        gzip -c ${fastq} > ${meta.id}.fastq.gz
        ln -s ${fastq} ${meta.id}.fastq
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(gzip --version 2>&1 | head -1)
    END_VERSIONS
    """
}
