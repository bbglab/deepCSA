include { SITESFROMPOSITIONS                                            } from '../../../modules/local/sitesfrompositions/main'
include { VCF_ANNOTATE_ENSEMBLVEP       as VCFANNOTATEPANEL             } from '../../../subworkflows/nf-core/vcf_annotate_ensemblvep_panel/main'
include { POSTPROCESS_VEP_ANNOTATION    as POSTPROCESSVEPPANEL          } from '../../../modules/local/process_annotation/main'

include { CREATECAPTUREDPANELS                                          } from '../../../modules/local/createpanels/captured/main'

include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSALL                  } from '../../../modules/local/createpanels/sample/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSPROTAFFECT           } from '../../../modules/local/createpanels/sample/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSNONPROTAFFECT        } from '../../../modules/local/createpanels/sample/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSEXONS                } from '../../../modules/local/createpanels/sample/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSINTRONS              } from '../../../modules/local/createpanels/sample/main'


include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSALL            } from '../../../modules/local/createpanels/consensus/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSPROTAFFECT     } from '../../../modules/local/createpanels/consensus/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSNONPROTAFFECT  } from '../../../modules/local/createpanels/consensus/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSEXONS          } from '../../../modules/local/createpanels/consensus/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSINTRONS        } from '../../../modules/local/createpanels/consensus/main'



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
    depths
    vep_cache
    vep_extra_files

    main:

    ch_versions = Channel.empty()

    // Create all possible sites and mutations per site of the captured panel
    SITESFROMPOSITIONS(depths)
    ch_versions = ch_versions.mix(SITESFROMPOSITIONS.out.versions)

    // Create a tuple for VEP annotation (mandatory)
    SITESFROMPOSITIONS.out.annotated_panel_reg.map{ it -> [[ id : "captured_panel"],  it[1]] }.set{ sites_annotation } // change names of vars
    // TEMPORARY: bgreference won't work in my cluster so to continue I give sites_annotation from params
    // Channel.of([1]).map{ it -> [[ id : "captured_panel"], params.sites_annotation] }.set{ sites_annotation }

    // Annotate all possible mutations in the captured panel
    // COMMENTED to avoid running VEP again and again during testing
    VCFANNOTATEPANEL(sites_annotation,
                    params.fasta,
                    params.vep_genome,
                    params.vep_species,
                    params.vep_cache_version,
                    vep_cache,
                    vep_extra_files)
    ch_versions = ch_versions.mix(VCFANNOTATEPANEL.out.versions)

    // Postprocess annotations to get one annotation per mutation
    POSTPROCESSVEPPANEL(VCFANNOTATEPANEL.out.tab)
    ch_versions = ch_versions.mix(POSTPROCESSVEPPANEL.out.versions)

    // Create captured-specific panels: all modalities
    CREATECAPTUREDPANELS(POSTPROCESSVEPPANEL.out.compact_panel_annotation)
    ch_versions = ch_versions.mix(CREATECAPTUREDPANELS.out.versions)

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


    // Create sample-specific panels: all modalities
    CREATESAMPLEPANELSALL(all_panel, depths, params.sample_panel_min_depth)
    CREATESAMPLEPANELSPROTAFFECT(prot_panel, depths, params.sample_panel_min_depth)
    CREATESAMPLEPANELSNONPROTAFFECT(nonprot_panel, depths, params.sample_panel_min_depth)
    CREATESAMPLEPANELSEXONS(exons_panel, depths, params.sample_panel_min_depth)
    CREATESAMPLEPANELSINTRONS(introns_panel, depths, params.sample_panel_min_depth)


    // Create consensus panel: all modalities
    CREATECONSENSUSPANELSALL(all_panel, depths, params.consensus_panel_min_depth)
    CREATECONSENSUSPANELSPROTAFFECT(prot_panel, depths, params.consensus_panel_min_depth)
    CREATECONSENSUSPANELSNONPROTAFFECT(nonprot_panel, depths, params.consensus_panel_min_depth)
    CREATECONSENSUSPANELSEXONS(exons_panel, depths, params.consensus_panel_min_depth)
    CREATECONSENSUSPANELSINTRONS(introns_panel, depths, params.consensus_panel_min_depth)


    emit:

    all_panel               = all_panel
    all_bed                 = all_bed
    prot_panel              = prot_panel
    prot_bed                = prot_bed
    nonprot_panel           = nonprot_panel
    nonprot_bed             = nonprot_bed
    exons_panel             = exons_panel
    exons_bed               = exons_bed
    introns_panel           = introns_panel
    introns_bed             = introns_bed

    all_consensus_panel     = CREATECONSENSUSPANELSALL.out.consensus_panel
    all_consensus_bed       = CREATECONSENSUSPANELSALL.out.consensus_panel_bed
    prot_consensus_panel    = CREATECONSENSUSPANELSPROTAFFECT.out.consensus_panel
    prot_consensus_bed      = CREATECONSENSUSPANELSPROTAFFECT.out.consensus_panel_bed
    nonprot_consensus_panel = CREATECONSENSUSPANELSNONPROTAFFECT.out.consensus_panel
    nonprot_consensus_bed   = CREATECONSENSUSPANELSNONPROTAFFECT.out.consensus_panel_bed
    exons_consensus_panel   = CREATECONSENSUSPANELSEXONS.out.consensus_panel
    exons_consensus_bed     = CREATECONSENSUSPANELSEXONS.out.consensus_panel_bed
    introns_consensus_panel = CREATECONSENSUSPANELSINTRONS.out.consensus_panel
    introns_consensus_bed   = CREATECONSENSUSPANELSINTRONS.out.consensus_panel_bed

    all_sample_panel        = restructureSamplePanel(CREATESAMPLEPANELSALL.out.sample_specific_panel.flatten())
    all_sample_bed          = restructureSamplePanel(CREATESAMPLEPANELSALL.out.sample_specific_panel_bed.flatten())
    prot_sample_panel       = restructureSamplePanel(CREATESAMPLEPANELSPROTAFFECT.out.sample_specific_panel.flatten())
    prot_sample_bed         = restructureSamplePanel(CREATESAMPLEPANELSPROTAFFECT.out.sample_specific_panel_bed.flatten())
    nonprot_sample_panel    = restructureSamplePanel(CREATESAMPLEPANELSNONPROTAFFECT.out.sample_specific_panel.flatten())
    nonprot_sample_bed      = restructureSamplePanel(CREATESAMPLEPANELSNONPROTAFFECT.out.sample_specific_panel_bed.flatten())
    exons_sample_panel      = restructureSamplePanel(CREATESAMPLEPANELSEXONS.out.sample_specific_panel.flatten())
    exons_sample_bed        = restructureSamplePanel(CREATESAMPLEPANELSEXONS.out.sample_specific_panel_bed.flatten())
    introns_sample_panel    = restructureSamplePanel(CREATESAMPLEPANELSINTRONS.out.sample_specific_panel.flatten())
    introns_sample_bed      = restructureSamplePanel(CREATESAMPLEPANELSINTRONS.out.sample_specific_panel_bed.flatten())

    versions = ch_versions                // channel: [ versions.yml ]

}
