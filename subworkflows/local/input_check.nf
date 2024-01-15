include { CHECK_SAMPLESHEET } from '../../modules/local/check_samplesheet'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    CHECK_SAMPLESHEET ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = CHECK_SAMPLESHEET.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = true

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read FastQ file does not exist!\n${row.fastq}"
    }
    fastq_meta = [ meta, [ file(row.fastq) ] ]

    return fastq_meta
}
