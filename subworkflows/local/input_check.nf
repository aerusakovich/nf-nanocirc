/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nanocirc input_check subworkflow
    Reads samplesheet CSV and emits [ meta, fastq ] tuples
    Expected columns: sample, fastq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INPUT_CHECK {

    take:
    samplesheet // file: path to samplesheet CSV

    main:
    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            // validate required columns
            if (!row.sample) exit 1, "Samplesheet missing 'sample' column: ${row}"
            if (!row.fastq)  exit 1, "Samplesheet missing 'fastq' column: ${row}"

            def meta  = [id: row.sample]
            def fastq = file(row.fastq, checkIfExists: true)
            [ meta, fastq ]
        }
        .set { fastq }

    emit:
    fastq // channel: [ val(meta), path(fastq) ]
}
