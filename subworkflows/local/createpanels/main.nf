include { SITESFROMPOSITIONS }                                     from '../../../modules/local/sitesfrompositions/main'
include { VCF_ANNOTATE_ALL              as VCFANNOTATEPANEL    }   from '../../../subworkflows/local/annotatepanel/main'
include { POSTPROCESS_VEP_ANNOTATION    as POSTPROCESSVEPPANEL }   from '../../../modules/local/process_annotation/main'
include { CREATECAPTUREDPANELS }                                   from '../../../modules/local/createcapturedpanels/main'
include { CREATESAMPLEPANELS }                                     from '../../../modules/local/createsamplepanels/main'


workflow CREATE_PANELS{

    take:
    depths
    vep_cache
    vep_extra_files

    main:

    ch_versions = Channel.empty()

    // Create all possible sites and mutations per site of the captured panel
    // SITESFROMPOSITIONS(depths)
    // ch_versions = ch_versions.mix(SITESFROMPOSITIONS.out.versions)

    // Create a tuple for VEP annotation (mandatory)
    // SITESFROMPOSITIONS.out.annotated_panel_reg.map{ it -> [[ id : "target_bed"],  it] }.set{ sites_annotation } //change names of vars
    // TEMPORARY: bgreference won't work in my cluster so to continue I give sites_annotation from params
    Channel.of([1]).map{ it -> [[ id : "captured_panel"], params.sites_annotation] }.set{ sites_annotation }

    // Annotate all possible mutations in the captured panel
    // COMMENTED to avoid running VEP again and again during testing
    // VCFANNOTATEPANEL(sites_annotation,
                    // params.fasta,
                    // params.vep_genome,
                    // params.vep_species,
                    // params.vep_cache_version,
                    // vep_cache,
                    // vep_extra_files)

    // Postprocess annotations to get one annotation per mutation
    // POSTPROCESSVEPPANEL(VCFANNOTATEPANEL.out.tab_ann)

    // Create captured-specific panels: all modalities
    // TEMPORARY: bgreference won't work in my cluster so to continue I give VCFANNOTATEPANEL.out.tab_ann from params
    // CREATECAPTUREDPANELS(POSTPROCESSVEPPANEL.out.compact_panel_annotation)
    Channel.of([1]).map{ it -> [[ id : "captured_panel"], params.compact_panel_annotation] }.set{ compact_panel_annotation }
    CREATECAPTUREDPANELS(compact_panel_annotation)

    // Create sample-specific panels: all modalities
    CREATECAPTUREDPANELS.out.collect().set{ captured_panels }
    CREATESAMPLEPANELS(captured_panels, depths, params.min_depth)

    // Create consensus panel: all modalities
    // CREATECONSENSUSPANEL()

    emit:
    // annotated_panel  = CREATECAPTUREDPANELS.out.captured_panel_protein_affecting   // channel: [ val(meta), file(depths) ]
    versions = ch_versions                // channel: [ versions.yml ]
}
