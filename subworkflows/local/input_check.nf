//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    require_bams

    main:
    SAMPLESHEET_CHECK ( samplesheet, require_bams)
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_input_channel(it, require_bams) }
        .set { sample_inputs }

    emit:
    sample_inputs                       // channel: [ val(meta), file(row.vcf), file(row.bam) ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_input_channel(LinkedHashMap row, required_bams=true) {
    // create meta map
    def meta = [:]
    meta.id       = row.sample
    meta.batch    = row.batch

    // add path(s) of the fastq file(s) to the meta map
    def vcf_bam_meta = []
    if (required_bams) {
        if (!file(row.vcf).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> VCF file does not exist!\n${row.vcf}"
        }
        if (!file(row.bam).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> BAM file does not exist!\n${row.bam}"
        }
        vcf_bam_meta = [ meta, file(row.vcf), file(row.bam) ]
    } else {
        if (!file(row.vcf).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> VCF file does not exist!\n${row.vcf}"
        }
        vcf_bam_meta = [ meta, file(row.vcf)]
    }

    

    return vcf_bam_meta
}
