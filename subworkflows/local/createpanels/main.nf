include { SITESFROMPOSITIONS                                        } from '../../../modules/local/sitesfrompositions/main'
include { VCF_ANNOTATE_ENSEMBLVEP       as VCFANNOTATEPANEL         } from '../../../subworkflows/nf-core/vcf_annotate_ensemblvep_panel/main'
include { POSTPROCESS_VEP_ANNOTATION    as POSTPROCESSVEPPANEL      } from '../../../modules/local/process_annotation/main'

include { CREATECAPTUREDPANELS                                      } from '../../../modules/local/createpanels/captured/main'

include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSALL              } from '../../../modules/local/createpanels/sample/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSPROTAFFECT       } from '../../../modules/local/createpanels/sample/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSNONPROTAFFECT    } from '../../../modules/local/createpanels/sample/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSEXONS            } from '../../../modules/local/createpanels/sample/main'
include { CREATESAMPLEPANELS as  CREATESAMPLEPANELSINTRONS          } from '../../../modules/local/createpanels/sample/main'


include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSALL        } from '../../../modules/local/createpanels/consensus/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSPROTAFFECT } from '../../../modules/local/createpanels/consensus/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSNONPROTAFFECT} from '../../../modules/local/createpanels/consensus/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSEXONS      } from '../../../modules/local/createpanels/consensus/main'
include { CREATECONSENSUSPANELS as  CREATECONSENSUSPANELSINTRONS    } from '../../../modules/local/createpanels/consensus/main'



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

    // Postprocess annotations to get one annotation per mutation
    POSTPROCESSVEPPANEL(VCFANNOTATEPANEL.out.tab)

    // Create captured-specific panels: all modalities
    CREATECAPTUREDPANELS(POSTPROCESSVEPPANEL.out.compact_panel_annotation)
    // Channel.of([1]).map{ it -> [[ id : "captured_panel"], params.compact_panel_annotation] }.set{ compact_panel_annotation }
    // CREATECAPTUREDPANELS(compact_panel_annotation)

    // Create sample-specific panels: all modalities
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
    all_panel               = CREATECAPTUREDPANELS.out.captured_panel_all
    all_bed                 = CREATECAPTUREDPANELS.out.captured_panel_all_bed
    prot_panel              = CREATECAPTUREDPANELS.out.captured_panel_protein_affecting
    prot_bed                = CREATECAPTUREDPANELS.out.captured_panel_protein_affecting_bed
    nonprot_panel           = CREATECAPTUREDPANELS.out.captured_panel_non_protein_affecting
    nonprot_bed             = CREATECAPTUREDPANELS.out.captured_panel_non_protein_affecting_bed
    exons_panel             = CREATECAPTUREDPANELS.out.captured_panel_exons_splice_sites
    exons_bed               = CREATECAPTUREDPANELS.out.captured_panel_exons_splice_sites_bed
    introns_panel           = CREATECAPTUREDPANELS.out.captured_panel_introns_intergenic
    introns_bed             = CREATECAPTUREDPANELS.out.captured_panel_introns_intergenic_bed


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


    all_sample_panel        = CREATESAMPLEPANELSALL.out.sample_specific_panel
    all_sample_bed          = CREATESAMPLEPANELSALL.out.sample_specific_panel_bed
    prot_sample_panel       = CREATESAMPLEPANELSPROTAFFECT.out.sample_specific_panel
    prot_sample_bed         = CREATESAMPLEPANELSPROTAFFECT.out.sample_specific_panel_bed
    nonprot_sample_panel    = CREATESAMPLEPANELSNONPROTAFFECT.out.sample_specific_panel
    nonprot_sample_bed      = CREATESAMPLEPANELSNONPROTAFFECT.out.sample_specific_panel_bed
    exons_sample_panel      = CREATESAMPLEPANELSEXONS.out.sample_specific_panel
    exons_sample_bed        = CREATESAMPLEPANELSEXONS.out.sample_specific_panel_bed
    introns_sample_panel    = CREATESAMPLEPANELSINTRONS.out.sample_specific_panel
    introns_sample_bed      = CREATESAMPLEPANELSINTRONS.out.sample_specific_panel_bed


    versions = ch_versions                // channel: [ versions.yml ]
}
