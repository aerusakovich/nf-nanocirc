process CIRCFL_SEQ {
    tag "$meta.id"
    label 'process_high'

    // --containall isolates the environment, --writable-tmpfs needed for temp files
    containerOptions = "--writable-tmpfs"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/anrusakovich/circfl-seq:latest' :
        'quay.io/anrusakovich/circfl-seq:latest' }"

    input:
    tuple val(meta), path(fastq)
    path  fasta
    path  gtf

    output:
    tuple val(meta), path("${meta.id}_circfl/circFL_final.bed"), emit: bed
    tuple val(meta), path("${meta.id}_circfl/"),                 emit: output_dir
    path  "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def out  = "${meta.id}_circfl"
    """
    export PYTHONNOUSERSITE=1

    # circfull requires a bgzip-compressed, tabix-indexed GTF sorted by position
    grep "^#" ${gtf} > anno_sorted.gtf
    grep -v "^#" ${gtf} | sort -k1,1 -k4,4n >> anno_sorted.gtf
    bgzip -c anno_sorted.gtf > anno.gtf.gz
    tabix -p gff anno.gtf.gz

    mkdir -p ${out}

    # Step 1: reference-guided detection
    circfull RG \\
        -f ${fastq} \\
        -g ${fasta} \\
        -a anno.gtf.gz \\
        -t ${task.cpus} \\
        -r \\
        -o ${out}

    # Step 2: de novo self-correction
    circfull DNSC \\
        -f ${out}/RG/circSeq.fa \\
        -t ${task.cpus} \\
        -o ${out}

    # Step 3: corrected reference-guided detection
    circfull cRG \\
        -f ${out}/DNSC \\
        -g ${fasta} \\
        -a anno.gtf.gz \\
        -t ${task.cpus} \\
        -o ${out}

    # Step 4: merge RG and cRG results
    circfull mRG \\
        -f ${fastq} \\
        -g ${fasta} \\
        -r ${out} \\
        -c ${out} \\
        -o ${out}

    # Step 5: annotate full-length circRNAs
    circfull anno \\
        -b ${out}/mRG/circFL_Normal_pass.bed \\
        -a anno.gtf.gz \\
        -o ${out}/annotated

    # Step 6: convert to BED12
    circfull FL2BED \\
        -c ${out}/mRG/circFL_Normal_pass.txt \\
        -o ${out}/circFL_final.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circfl_seq: \$(circfull --version 2>&1 | head -1)
    END_VERSIONS
    """
}
