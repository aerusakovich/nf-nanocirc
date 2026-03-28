process PREPARE_GENOME {
    tag "prepare_genome"
    label 'process_low'

    // Copy FASTA and build BWA + samtools fai indices once.
    // storeDir caches results — skipped if files already exist.
    storeDir "${params.genome_index_dir ?: "${params.outdir}/genome_index"}"

    // bwakit contains both bwa and samtools
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwakit:0.7.15--he513fc3_2' :
        'quay.io/biocontainers/bwakit:0.7.15--he513fc3_2' }"

    input:
    path fasta

    output:
    path "genome.fa",    emit: fasta   // real file — bwapy finds index alongside it
    path "genome.fa.*",  emit: index   // bwa + samtools fai indices
    path "versions.yml", emit: versions

    script:
    """
    # Resolve symlink so bwapy finds index next to the real genome.fa
    cp -L ${fasta} genome.fa

    # BWA index (required by CIRI-long via bwapy)
    bwa index genome.fa

    # samtools fai index (required by CIRI-long and other tools)
    samtools faidx genome.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep Version | sed 's/Version: //')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
