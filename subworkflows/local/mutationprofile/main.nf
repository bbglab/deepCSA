include { COMPUTEDEPTHS } from '../../modules/local/computedepths/main'


workflow MUTATIONAL_PROFILE{
    take:
    // inputs
    bam_list

    main:
    // actual code
    ch_versions = Channel.empty()
    COMPUTEDEPTHS(bam_list)

    // PROCESSDEPTHSTABLE()

    // PLOTDEPTHS()

    emit:
    depths   = COMPUTEDEPTHS.out.depths   // channel: [ val(meta), file(depths) ]
    versions = ch_versions                // channel: [ versions.yml ]
}
