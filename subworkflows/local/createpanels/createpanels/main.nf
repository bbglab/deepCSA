include { CREATECAPTUREDPANELS                                          } from '../../../../modules/local/createpanels/captured/main'

include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSALL                  } from '../../../../modules/local/createpanels/sample/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSPROTAFFECT           } from '../../../../modules/local/createpanels/sample/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSNONPROTAFFECT        } from '../../../../modules/local/createpanels/sample/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSEXONS                } from '../../../../modules/local/createpanels/sample/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSINTRONS              } from '../../../../modules/local/createpanels/sample/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSSYNONYMOUS           } from '../../../../modules/local/createpanels/sample/main'


include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSALL            } from '../../../../modules/local/createpanels/consensus/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSPROTAFFECT     } from '../../../../modules/local/createpanels/consensus/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSNONPROTAFFECT  } from '../../../../modules/local/createpanels/consensus/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSEXONS          } from '../../../../modules/local/createpanels/consensus/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSINTRONS        } from '../../../../modules/local/createpanels/consensus/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSSYNONYMOUS     } from '../../../../modules/local/createpanels/consensus/main'



def restructurePanel(panel) {
        // return panel.map{ it -> [[id: it[0].name.tokenize('.')[0..1].join('.')], it[1]] }
        // return panel.map { it -> [[id: it.name.tokenize('.')[0..1].join('.')], it] }
        return panel.map { it -> [[id: it.name.tokenize('.')[1]], it] }
    }


def restructureSamplePanel(panel) {
        // return panel.map{ it -> [[id: it[0].name.tokenize('.')[0..1].join('.')], it[1]] }
        // return panel.map { it -> [[id: it.name.tokenize('.')[0..1].join('.')], it] }
        return panel.map { it -> [[id: it.name.tokenize('.')[0]], it] }
    }



workflow CREATE_PANELS {

    take:
    complete_annotated_panel
    depths
    // flagged_regions_bed  // BED file with regions to exclude -> from MUTATION_PREPROCESSING

    main:

    // // Filter out flagged regions from the annotated panel before creating panels
    // // This happens via bedtools subtract or similar filtering in the annotation
    // if (flagged_regions_bed) {
    //     // TODO: Add a process here to subtract flagged_regions_bed from complete_annotated_panel
    //     // For now, just use the panel as-is (implement filtering process separately)
    //     filtered_panel = complete_annotated_panel
    // } else {
    //     filtered_panel = complete_annotated_panel
    // }

    // Create captured-specific panels: all modalities
    CREATECAPTUREDPANELS(complete_annotated_panel)

    restructurePanel(CREATECAPTUREDPANELS.out.captured_panel_all).set{all_panel}
    restructurePanel(CREATECAPTUREDPANELS.out.captured_panel_all_bed).set{all_bed}
    restructurePanel(CREATECAPTUREDPANELS.out.captured_panel_protein_affecting).set{prot_panel}
    restructurePanel(CREATECAPTUREDPANELS.out.captured_panel_protein_affecting_bed).set{prot_bed}
    restructurePanel(CREATECAPTUREDPANELS.out.captured_panel_non_protein_affecting).set{nonprot_panel}
    restructurePanel(CREATECAPTUREDPANELS.out.captured_panel_non_protein_affecting_bed).set{nonprot_bed}
    restructurePanel(CREATECAPTUREDPANELS.out.captured_panel_exons_splice_sites).set{exons_panel}
    restructurePanel(CREATECAPTUREDPANELS.out.captured_panel_exons_splice_sites_bed).set{exons_bed}
    restructurePanel(CREATECAPTUREDPANELS.out.captured_panel_introns_intergenic).set{introns_panel}
    restructurePanel(CREATECAPTUREDPANELS.out.captured_panel_introns_intergenic_bed).set{introns_bed}
    restructurePanel(CREATECAPTUREDPANELS.out.captured_panel_synonymous).set{synonymous_panel}
    restructurePanel(CREATECAPTUREDPANELS.out.captured_panel_synonymous_bed).set{synonymous_bed}

    if (params.create_sample_panels){
        // Create sample-specific panels: all modalities
        CREATESAMPLEPANELSALL(all_panel, depths, params.sample_panel_min_depth)
        CREATESAMPLEPANELSPROTAFFECT(prot_panel, depths, params.sample_panel_min_depth)
        CREATESAMPLEPANELSNONPROTAFFECT(nonprot_panel, depths, params.sample_panel_min_depth)
        CREATESAMPLEPANELSEXONS(exons_panel, depths, params.sample_panel_min_depth)
        CREATESAMPLEPANELSINTRONS(introns_panel, depths, params.sample_panel_min_depth)
        CREATESAMPLEPANELSSYNONYMOUS(synonymous_panel, depths, params.sample_panel_min_depth)
    }


    // Create consensus panel: all modalities
    CREATECONSENSUSPANELSALL(all_panel, depths)
    CREATECONSENSUSPANELSPROTAFFECT(prot_panel, depths)
    CREATECONSENSUSPANELSNONPROTAFFECT(nonprot_panel, depths)
    CREATECONSENSUSPANELSEXONS(exons_panel, depths)
    CREATECONSENSUSPANELSINTRONS(introns_panel, depths)
    CREATECONSENSUSPANELSSYNONYMOUS(synonymous_panel, depths)

    emit:
    // full_panel_annotated        = VCFANNOTATEPANEL.out.tab
    all_panel                   = all_panel.first()
    all_bed                     = all_bed.first()
    prot_panel                  = prot_panel.first()
    prot_bed                    = prot_bed.first()
    nonprot_panel               = nonprot_panel.first()
    nonprot_bed                 = nonprot_bed.first()
    exons_panel                 = exons_panel.first()
    exons_bed                   = exons_bed.first()
    introns_panel               = introns_panel.first()
    introns_bed                 = introns_bed.first()
    synonymous_panel            = synonymous_panel.first()
    synonymous_bed              = synonymous_bed.first()


    all_consensus_panel         = CREATECONSENSUSPANELSALL.out.consensus_panel.first()
    all_consensus_bed           = CREATECONSENSUSPANELSALL.out.consensus_panel_bed.first()
    prot_consensus_panel        = CREATECONSENSUSPANELSPROTAFFECT.out.consensus_panel.first()
    prot_consensus_bed          = CREATECONSENSUSPANELSPROTAFFECT.out.consensus_panel_bed.first()
    nonprot_consensus_panel     = CREATECONSENSUSPANELSNONPROTAFFECT.out.consensus_panel.first()
    nonprot_consensus_bed       = CREATECONSENSUSPANELSNONPROTAFFECT.out.consensus_panel_bed.first()
    exons_consensus_panel       = CREATECONSENSUSPANELSEXONS.out.consensus_panel.first()
    exons_consensus_bed         = CREATECONSENSUSPANELSEXONS.out.consensus_panel_bed.first()
    introns_consensus_panel     = CREATECONSENSUSPANELSINTRONS.out.consensus_panel.first()
    introns_consensus_bed       = CREATECONSENSUSPANELSINTRONS.out.consensus_panel_bed.first()
    synonymous_consensus_panel  = CREATECONSENSUSPANELSSYNONYMOUS.out.consensus_panel.first()
    synonymous_consensus_bed    = CREATECONSENSUSPANELSSYNONYMOUS.out.consensus_panel_bed.first()

    // all_sample_panel        = restructureSamplePanel(CREATESAMPLEPANELSALL.out.sample_specific_panel.flatten())
    // all_sample_bed          = restructureSamplePanel(CREATESAMPLEPANELSALL.out.sample_specific_panel_bed.flatten())
    // prot_sample_panel       = restructureSamplePanel(CREATESAMPLEPANELSPROTAFFECT.out.sample_specific_panel.flatten())
    // prot_sample_bed         = restructureSamplePanel(CREATESAMPLEPANELSPROTAFFECT.out.sample_specific_panel_bed.flatten())
    // nonprot_sample_panel    = restructureSamplePanel(CREATESAMPLEPANELSNONPROTAFFECT.out.sample_specific_panel.flatten())
    // nonprot_sample_bed      = restructureSamplePanel(CREATESAMPLEPANELSNONPROTAFFECT.out.sample_specific_panel_bed.flatten())
    // exons_sample_panel      = restructureSamplePanel(CREATESAMPLEPANELSEXONS.out.sample_specific_panel.flatten())
    // exons_sample_bed        = restructureSamplePanel(CREATESAMPLEPANELSEXONS.out.sample_specific_panel_bed.flatten())
    // introns_sample_panel    = restructureSamplePanel(CREATESAMPLEPANELSINTRONS.out.sample_specific_panel.flatten())
    // introns_sample_bed      = restructureSamplePanel(CREATESAMPLEPANELSINTRONS.out.sample_specific_panel_bed.flatten())

}
