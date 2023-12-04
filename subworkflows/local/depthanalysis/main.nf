include { COMPUTEDEPTHS } from '../../modules/local/computedepths/main'


workflow DEPTH_ANALYSIS{
    take:
    bam_list

    main:
    COMPUTEDEPTHS(bam)

    // PROCESSDEPTHSTABLE()

    // PLOTDEPTHS()

    emit:
    COMPUTEDEPTHS.out.depths                  // channel: [ val(meta), file(row.vcf), file(row.bam) ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}