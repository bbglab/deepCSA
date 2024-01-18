include { SITESFROMPOSITIONS }                                           from '../../../modules/local/sitesfrompositions/main'
include { VCF_ANNOTATE_ALL              as VCFANNOTATEPANEL    }         from '../../../subworkflows/local/annotatepanel/main'
include { POSTPROCESS_VEP_ANNOTATION    as POSTPROCESSVEPPANEL }         from '../../../modules/local/process_annotation/main'

include { CREATECAPTUREDPANELS }                                         from '../../../modules/local/createcapturedpanels/main'

include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSALL}                  from '../../../modules/local/createsamplepanels/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSPROTAFFECT}           from '../../../modules/local/createsamplepanels/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSNONPROTAFFECT}        from '../../../modules/local/createsamplepanels/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSEXONS}                from '../../../modules/local/createsamplepanels/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSINTRONS}              from '../../../modules/local/createsamplepanels/main'


include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSALL}              from '../../../modules/local/createconsensuspanels/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSPROTAFFECT}       from '../../../modules/local/createconsensuspanels/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSNONPROTAFFECT}    from '../../../modules/local/createconsensuspanels/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSEXONS}            from '../../../modules/local/createconsensuspanels/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSINTRONS}          from '../../../modules/local/createconsensuspanels/main'



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
    // CREATECAPTUREDPANELS.out.captured_panel_protein_affecting.concat(
	// 										CREATECAPTUREDPANELS.out.captured_panel_non_protein_affecting,
	// 										CREATECAPTUREDPANELS.out.captured_panel_exons_splice_sites,
	// 										CREATECAPTUREDPANELS.out.captured_panel_introns_intergenic )
	// 										.set{ captured_panels }
    CREATESAMPLEPANELSALL(CREATECAPTUREDPANELS.out.captured_panel_all, depths, params.min_depth)
    CREATESAMPLEPANELSPROTAFFECT(CREATECAPTUREDPANELS.out.captured_panel_protein_affecting, depths, params.min_depth)
    CREATESAMPLEPANELSNONPROTAFFECT(CREATECAPTUREDPANELS.out.captured_panel_non_protein_affecting, depths, params.min_depth)
    CREATESAMPLEPANELSEXONS(CREATECAPTUREDPANELS.out.captured_panel_exons_splice_sites, depths, params.min_depth)
    CREATESAMPLEPANELSINTRONS(CREATECAPTUREDPANELS.out.captured_panel_introns_intergenic, depths, params.min_depth)

    // Create consensus panel: all modalities
    CREATECONSENSUSPANELSALL(CREATECAPTUREDPANELS.out.captured_panel_all, depths, params.consensus_min_depth)
    CREATECONSENSUSPANELSPROTAFFECT(CREATECAPTUREDPANELS.out.captured_panel_protein_affecting, depths, params.consensus_min_depth)
    CREATECONSENSUSPANELSNONPROTAFFECT(CREATECAPTUREDPANELS.out.captured_panel_non_protein_affecting, depths, params.consensus_min_depth)
    CREATECONSENSUSPANELSEXONS(CREATECAPTUREDPANELS.out.captured_panel_exons_splice_sites, depths, params.consensus_min_depth)
    CREATECONSENSUSPANELSINTRONS(CREATECAPTUREDPANELS.out.captured_panel_introns_intergenic, depths, params.consensus_min_depth)

    emit:
    // annotated_panel  = CREATECAPTUREDPANELS.out.captured_panel_protein_affecting   // channel: [ val(meta), file(depths) ]
    versions = ch_versions                // channel: [ versions.yml ]
}
