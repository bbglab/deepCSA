include { CREATECAPTUREDPANELS                                          } from '../../../../modules/local/createpanels/captured/main'

include { FILTER_CAPTURED_PANEL as  FILTERCAPTUREDPANEL                 } from '../../../../modules/local/filter_captured_panel/main'

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

workflow CREATE_PANELS {

    take:
    complete_annotated_panel
    depths
    flagged_bed  // BED file with regions to exclude -> from MUT_PREPROCESSING

    main:

    // Filter out flagged regions from the annotated panel before creating panels - bedtools intersect
    FILTERCAPTUREDPANEL(complete_annotated_panel, flagged_bed)
    filtered_panel = FILTERCAPTUREDPANEL.out.filtered_annotated_panel

    // Create captured-specific panels: all modalities
    CREATECAPTUREDPANELS(filtered_panel)

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

    // Create consensus panel: all modalities
    CREATECONSENSUSPANELSALL(all_panel, depths)
    CREATECONSENSUSPANELSPROTAFFECT(prot_panel, depths)
    CREATECONSENSUSPANELSNONPROTAFFECT(nonprot_panel, depths)
    CREATECONSENSUSPANELSEXONS(exons_panel, depths)
    CREATECONSENSUSPANELSINTRONS(introns_panel, depths)
    CREATECONSENSUSPANELSSYNONYMOUS(synonymous_panel, depths)

    emit:
    removed_sites               = FILTERCAPTUREDPANEL.out.removed_sites

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
}
