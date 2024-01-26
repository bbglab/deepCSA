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
    // TODO confirm if this is true
    // here every process generates a multiple tsvs and beds, one for each sample
    CREATESAMPLEPANELSALL(CREATECAPTUREDPANELS.out.captured_panel_all, depths, params.min_depth)
    CREATESAMPLEPANELSPROTAFFECT(CREATECAPTUREDPANELS.out.captured_panel_protein_affecting, depths, params.min_depth)
    CREATESAMPLEPANELSNONPROTAFFECT(CREATECAPTUREDPANELS.out.captured_panel_non_protein_affecting, depths, params.min_depth)
    CREATESAMPLEPANELSEXONS(CREATECAPTUREDPANELS.out.captured_panel_exons_splice_sites, depths, params.min_depth)
    CREATESAMPLEPANELSINTRONS(CREATECAPTUREDPANELS.out.captured_panel_introns_intergenic, depths, params.min_depth)

    // Create consensus panel: all modalities
    // TODO confirm if this is true
    // here every process generates a single tsv and bed file
    CREATECONSENSUSPANELSALL(CREATECAPTUREDPANELS.out.captured_panel_all, depths, params.consensus_min_depth)
    CREATECONSENSUSPANELSPROTAFFECT(CREATECAPTUREDPANELS.out.captured_panel_protein_affecting, depths, params.consensus_min_depth)
    CREATECONSENSUSPANELSNONPROTAFFECT(CREATECAPTUREDPANELS.out.captured_panel_non_protein_affecting, depths, params.consensus_min_depth)
    CREATECONSENSUSPANELSEXONS(CREATECAPTUREDPANELS.out.captured_panel_exons_splice_sites, depths, params.consensus_min_depth)
    CREATECONSENSUSPANELSINTRONS(CREATECAPTUREDPANELS.out.captured_panel_introns_intergenic, depths, params.consensus_min_depth)

    emit:
    annotated_panel         = CREATECAPTUREDPANELS.out.captured_panel_protein_affecting   // channel: [ val(meta), file(depths) ]
    exons_consensus_panel   = CREATECONSENSUSPANELSEXONS.out.consensus_panel
    exons_consensus_bed     = CREATECONSENSUSPANELSEXONS.out.consensus_panel_bed
    captured_panel          = CREATECAPTUREDPANELS.out.captured_panel_all
    captured_bed            = CREATECAPTUREDPANELS.out.captured_panel_all_bed
    introns_consensus_panel = CREATECONSENSUSPANELSINTRONS.out.consensus_panel
    introns_consensus_bed   = CREATECONSENSUSPANELSINTRONS.out.consensus_panel_bed
    // captured_bed            = CREATESAMPLEPANELSNONPROTAFFECT.out.captured_panel_all_bed
    versions = ch_versions                // channel: [ versions.yml ]
}
