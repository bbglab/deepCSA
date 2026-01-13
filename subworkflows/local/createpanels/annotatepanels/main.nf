include { SITESFROMPOSITIONS                                            } from '../../../../modules/local/sitesfrompositions/main'
include { VCF_ANNOTATE_ENSEMBLVEP       as VCFANNOTATEPANEL             } from '../../../../subworkflows/nf-core/vcf_annotate_ensemblvep_panel/main'

include { POSTPROCESS_VEP_ANNOTATION    as POSTPROCESSVEPPANEL          } from '../../../../modules/local/process_annotation/panel/main'

include { CUSTOM_ANNOTATION_PROCESSING  as CUSTOMPROCESSING             } from '../../../../modules/local/process_annotation/panelcustom/main'
include { CUSTOM_ANNOTATION_PROCESSING  as CUSTOMPROCESSINGRICH         } from '../../../../modules/local/process_annotation/panelcustom/main'

include { DOMAIN_ANNOTATION             as DOMAINANNOTATION             } from '../../../../modules/local/process_annotation/domain/main'

include { CREATECAPTUREDPANELS                                          } from '../../../../modules/local/createpanels/captured/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSALL                  } from '../../../../modules/local/createpanels/sample/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSALL            } from '../../../../modules/local/createpanels/consensus/main'

def restructurePanel(panel) {
        // return panel.map{ it -> [[id: it[0].name.tokenize('.')[0..1].join('.')], it[1]] }
        // return panel.map { it -> [[id: it.name.tokenize('.')[0..1].join('.')], it] }
        return panel.map { it -> [[id: it.name.tokenize('.')[1]], it] }
    }


workflow ANNOTATE_PANELS {

    take:
    depths

    main:

    // Create all possible sites and mutations per site of the captured panel
    SITESFROMPOSITIONS(depths)

    // Create a tuple for VEP annotation (mandatory)
    SITESFROMPOSITIONS.out.annotated_panel_reg.map{ it -> [[ id : "captured_panel"],  it[1]] }.set{ sites_annotation }

    // Annotate all possible mutations in the captured panel
    VCFANNOTATEPANEL(sites_annotation,
                    params.fasta,
                    params.vep_genome,
                    params.vep_species,
                    params.vep_cache_version,
                    params.vep_cache,
                    [])

    // Postprocess annotations to get one annotation per mutation
    POSTPROCESSVEPPANEL(VCFANNOTATEPANEL.out.tab)

    if (params.customize_annotation) {
        custom_annotation_tsv = file(params.custom_annotation_tsv)

        // Update specific regions based on user preferences
        CUSTOMPROCESSING(POSTPROCESSVEPPANEL.out.compact_panel_annotation, custom_annotation_tsv)
        complete_annotated_panel = CUSTOMPROCESSING.out.custom_panel_annotation

        CUSTOMPROCESSINGRICH(POSTPROCESSVEPPANEL.out.rich_panel_annotation, custom_annotation_tsv)
        rich_annotated = CUSTOMPROCESSINGRICH.out.custom_panel_annotation

        added_regions = CUSTOMPROCESSINGRICH.out.added_regions

    } else {
        complete_annotated_panel = POSTPROCESSVEPPANEL.out.compact_panel_annotation
        rich_annotated = POSTPROCESSVEPPANEL.out.rich_panel_annotation
        added_regions = Channel.empty()
    }

    domains = file(params.domains_file, checkIfExists: true)
    DOMAINANNOTATION(rich_annotated, domains)

    // Create MUTPREPROCESSING-specific panels
    CREATECAPTUREDPANELS(complete_annotated_panel)

    // All consensus bed
    restructurePanel(CREATECAPTUREDPANELS.out.captured_panel_all).set{all_panel}
    if (params.create_sample_panels){
        CREATESAMPLEPANELSALL(all_panel, depths, params.sample_panel_min_depth)
    }
    CREATECONSENSUSPANELSALL(all_panel, depths)
    
    // Exons bed
    restructurePanel(CREATECAPTUREDPANELS.out.captured_panel_exons_splice_sites_bed).set{exons_bed_initial}
    
    emit:
    full_panel_annotated           = VCFANNOTATEPANEL.out.tab
    complete_annotated_panel       = complete_annotated_panel

    panel_annotated_rich           = rich_annotated
    added_custom_regions           = added_regions
    domains_panel_bed              = DOMAINANNOTATION.out.domains_bed.first()
    domains_in_panel               = DOMAINANNOTATION.out.domains_tsv.first()

    postprocessed_panel            = POSTPROCESSVEPPANEL.out.compact_panel_annotation.first()
    postprocessed_panel_rich       = POSTPROCESSVEPPANEL.out.rich_panel_annotation.first()

    all_consensus_bed_initial    = CREATECONSENSUSPANELSALL.out.consensus_panel_bed.first()    
    exons_bed_initial              = exons_bed_initial.first()
}
