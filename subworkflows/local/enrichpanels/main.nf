include { DNA_2_PROTEIN_MAPPING     as DNA2PROTEINMAPPING       }   from '../../../modules/local/dna2protein/main'

include { EXPAND_REGIONS            as EXPANDREGIONSALL         }   from '../../../modules/local/expand_regions/main'
include { EXPAND_REGIONS            as EXPANDREGIONSNONPROT     }   from '../../../modules/local/expand_regions/main'
include { EXPAND_REGIONS            as EXPANDREGIONSPROT        }   from '../../../modules/local/expand_regions/main'
include { EXPAND_REGIONS            as EXPANDREGIONSSYNONYMOUS  }   from '../../../modules/local/expand_regions/main'
include { EXPAND_REGIONS            as EXPANDREGIONSEXONS       }   from '../../../modules/local/expand_regions/main'



workflow ENRICHPANELS {

    take:
    mutations
    all_samples_depths

    all_consensus_panel
    nonprot_consensus_panel
    prot_consensus_panel
    synonymous_consensus_panel
    exons_consensus_panel

    domains_file

    exons_bedfile

    main:

    gff3_channel = params.gff3_file ? file(params.gff3_file, checkIfExists: true) : file(params.input, checkIfExists: true)
    DNA2PROTEINMAPPING(mutations, exons_consensus_panel, all_samples_depths, gff3_channel)

    // Create a channel for the domains file if autodomains is true
    domains_ch = params.autodomains ? domains_file : []  // .map{ it -> it[1]} : []
    exons_ch = params.autoexons ? DNA2PROTEINMAPPING.out.panel_exons_bed.map{ it -> it[1]}.ifEmpty([])  : []

    // Create a channel for the subgenic bedfile if provided
    subgenic_ch = params.subgenic_bedfile ? file(params.subgenic_bedfile, checkIfExists: true) : []

    if (params.create_subgenic_regions){
        EXPANDREGIONSALL(all_consensus_panel, domains_ch, exons_ch, subgenic_ch)
        all_expanded_panel = EXPANDREGIONSALL.out.panel_increased.first()

        EXPANDREGIONSNONPROT(nonprot_consensus_panel, domains_ch, exons_ch, subgenic_ch)
        nonprot_expanded_panel = EXPANDREGIONSNONPROT.out.panel_increased.first()

        EXPANDREGIONSPROT(prot_consensus_panel, domains_ch, exons_ch, subgenic_ch)
        prot_expanded_panel = EXPANDREGIONSPROT.out.panel_increased.first()

        EXPANDREGIONSSYNONYMOUS(synonymous_consensus_panel, domains_ch, exons_ch, subgenic_ch)
        synonymous_expanded_panel = EXPANDREGIONSSYNONYMOUS.out.panel_increased.first()

        EXPANDREGIONSEXONS(exons_consensus_panel, domains_ch, exons_ch, subgenic_ch)
        exons_expanded_panel = EXPANDREGIONSEXONS.out.panel_increased.first()

        // all_json_subgenic = EXPANDREGIONSALL.out.new_regions_json
        // nonprot_json_subgenic = EXPANDREGIONSNONPROT.out.new_regions_json
        // prot_json_subgenic = EXPANDREGIONSPROT.out.new_regions_json
        // synonymous_json_subgenic = EXPANDREGIONSSYNONYMOUS.out.new_regions_json
        exons_json_subgenic = EXPANDREGIONSEXONS.out.new_regions_json.first()

    } else {
        all_expanded_panel          = all_consensus_panel
        nonprot_expanded_panel      = nonprot_consensus_panel
        prot_expanded_panel         = prot_consensus_panel
        synonymous_expanded_panel   = synonymous_consensus_panel
        exons_expanded_panel        = exons_consensus_panel

        exons_json_subgenic         = exons_bedfile
    }

    emit:

    all_consensus_expanded_panel        = all_expanded_panel
    nonprot_consensus_expanded_panel    = nonprot_expanded_panel
    prot_consensus_expanded_panel       = prot_expanded_panel
    synonymous_consensus_expanded_panel = synonymous_expanded_panel
    exons_consensus_expanded_panel      = exons_expanded_panel

    exons_json_subgenic                 = exons_json_subgenic

    dna2protein_mapping_depth_exons     = DNA2PROTEINMAPPING.out.depths_exons_positions.first()
    dna2protein_mapping_panel_exons     = DNA2PROTEINMAPPING.out.panel_exons_bed.first()
}
