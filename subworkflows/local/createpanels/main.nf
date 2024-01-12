include { SITESFROMPOSITIONS } from '../../../modules/local/sitesfrompositions/main'
// include { VCF_ANNOTATE_ENSEMBLVEP   as VCFANNOTATEPANEL  } from '../../nf-core/vcf_annotate_ensemblvep/main'
include { VCF_ANNOTATE_ALL   as VCFANNOTATEPANEL    } from '../../../subworkflows/local/annotatepanel/main'


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
    Channel.of([1]).map{ it -> [[ id : "target_bed"], params.sites_annotation] }.set{ sites_annotation }

    // Annotate all possible mutations in the captured panel
    VCFANNOTATEPANEL(sites_annotation,
                                    params.fasta,
                                    params.vep_genome,
                                    params.vep_species,
                                    params.vep_cache_version,
                                    vep_cache,
                                    vep_extra_files)

    // // Postprocess annotations to get one annotation per mutation
    // POSTPROCESSVEPPANEL(VCFANNOTATEPANEL.out.tab_ann)

    // // Create sample-specific panels: all modalities
    // CREATESAMPLEPANEL()

    // // Create consensus panel: all modalities
    // CREATECONSENSUSPANEL()

    emit:
    annotated_panel  = VCFANNOTATEPANEL.out.tab_ann   // channel: [ val(meta), file(depths) ]
    versions = ch_versions                // channel: [ versions.yml ]
}
