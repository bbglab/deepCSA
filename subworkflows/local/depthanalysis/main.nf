include { COMPUTEDEPTHS } from '../../../modules/local/computedepths/main'


workflow DEPTH_ANALYSIS{
    take:
    // inputs
    bam_list

    main:
    // actual code

    ch_versions = Channel.empty()

    // Join all annotated samples and put them in a channel to be summarized together
    bam_list.map{ it -> it[1] }.collect().map{ it -> [[ id:"all_samples" ], it]}.set{ combined_bams }

    COMPUTEDEPTHS(combined_bams)

    // PROCESSDEPTHSTABLE()

    // PLOTDEPTHS()

    emit:
    depths   = COMPUTEDEPTHS.out.depths   // channel: [ val(meta), file(depths) ]
    versions = ch_versions                // channel: [ versions.yml ]
}
