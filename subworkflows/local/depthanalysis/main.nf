include { COMPUTEDEPTHS } from '../../../modules/local/computedepths/main'


workflow DEPTH_ANALYSIS{

    take:
    input_files
    custom_bed

    main:

    if (params.use_custom_depths) {
        output_depths = channel.fromPath( params.custom_depths_table, checkIfExists: true).map{ path -> [ [id: "all_samples"], path ] }.first()
    } else {
        input_files
        .map{ it -> [it[0], it[2]]}
        .set{ bam_list }

        // Join all samples and put them in a channel to be summarized together
        bam_list.map{ it -> it[1] }.collect().map{ it -> [[ id:"all_samples" ], it]}.set{ combined_bams }

        // Create a table with the depth per positions across samples
        COMPUTEDEPTHS(combined_bams, custom_bed)
        output_depths = COMPUTEDEPTHS.out.depths
    }

    emit:
    depths   = output_depths
}
