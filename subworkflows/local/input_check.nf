/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nanocirc input_check subworkflow
    Reads samplesheet CSV and emits [ meta, fastq ] tuples
    Expected columns: sample, fastq
    Input files must be gzipped (.fastq.gz or .fq.gz)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INPUT_CHECK {

    take:
    samplesheet

    main:
    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            if (!row.sample) exit 1, "Samplesheet missing 'sample' column: ${row}"
            if (!row.fastq)  exit 1, "Samplesheet missing 'fastq' column: ${row}"

            def meta  = [id: row.sample]
            def fastq = file(row.fastq, checkIfExists: true)

            if (!fastq.name.endsWith('.gz')) {
                exit 1, "Sample '${row.sample}': '${fastq.name}' must be gzipped (.fastq.gz or .fq.gz). Please compress your input files."
            }

            [ meta, fastq ]
        }
        .set { fastq }

    emit:
    fastq
}
