include { SITESFROMPOSITIONS } from '../../../modules/local/sitesfrompositions/main'


workflow CREATE_PANELS{

    take:
    depths

    main:

    ch_versions = Channel.empty()

    // Create all possible sites and mutations per site of the captured panel
    SITESFROMPOSITIONS(depths)
    ch_versions = ch_versions.mix(SITESFROMPOSITIONS.out.versions)

    // // Create a tuple for VEP annotation (mandatory)
    // SITESFROMPOSITIONS.out.annotated_panel_reg.map{ it -> [[ id : "target_bed"],  it] }.set{ sites_annotation } //change names of vars

    // // Annotate all possible mutations in the captured panel
    // VCFANNOTATEPANEL(sites_annotation,
    //                                 ch_ref_fasta,
    //                                 params.vep_genome,
    //                                 params.vep_species,
    //                                 params.vep_cache_version,
    //                                 vep_cache,
    //                                 vep_extra_files)

    // // Postprocess annotations to get one annotation per mutation
    // POSTPROCESSVEPPANEL(VCFANNOTATEPANEL.out.tab_ann)

    // // Create sample-specific panels: all modalities
    // CREATESAMPLEPANEL()

    // // Create consensus panel: all modalities
    // CREATECONSENSUSPANEL()

    emit:
    annotated_panel_reg   = SITESFROMPOSITIONS.out.annotated_panel_reg   // channel: [ val(meta), file(depths) ]
    versions = ch_versions                // channel: [ versions.yml ]
}
