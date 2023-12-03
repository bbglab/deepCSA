//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_input_channel(it) }
        .set { mutations }

    emit:
    mutations                                 // channel: [ val(meta), file(row.vcf), file(row.bam) ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_input_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id       = row.sample
    meta.batch    = row.batch
    // meta.vcf_only = row.vcf_only.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def vcf_bam_meta = []
    if (!file(row.vcf).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> VCF file does not exist!\n${row.vcf}"
    }
    if (!file(row.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> BAM file does not exist!\n${row.vcf}"
    }
    
    vcf_bam_meta = [ meta, file(row.vcf), file(row.bam) ]

    return vcf_bam_meta
}
