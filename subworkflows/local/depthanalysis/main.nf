include { COMPUTEDEPTHS } from '../../../modules/local/computedepths/main'


workflow DEPTH_ANALYSIS {
    take:
    bam_list
    custom_bed

    main:


    // Join all samples and put them in a channel to be summarized together
    bam_list.map { it -> it[1] }.collect().map { it -> [[id: "all_samples"], it] }.set { combined_bams }

    // Create a table with the depth per positions across samples
    COMPUTEDEPTHS(combined_bams, custom_bed)

    emit:
    depths = COMPUTEDEPTHS.out.depths // channel: [ val(meta), file(depths) ]
}
