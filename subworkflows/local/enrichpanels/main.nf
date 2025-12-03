include { DNA_2_PROTEIN_MAPPING     as DNA2PROTEINMAPPING               }   from '../../../modules/local/dna2protein/main'

include { EXPAND_REGIONS            as EXPANDREGIONSALL                    }   from '../../../modules/local/expand_regions/main'
include { EXPAND_REGIONS            as EXPANDREGIONSNONPROT                    }   from '../../../modules/local/expand_regions/main'
include { EXPAND_REGIONS            as EXPANDREGIONSPROT                    }   from '../../../modules/local/expand_regions/main'
include { EXPAND_REGIONS            as EXPANDREGIONSSYNONYMOUS                    }   from '../../../modules/local/expand_regions/main'
include { EXPAND_REGIONS            as EXPANDREGIONSEXONS                    }   from '../../../modules/local/expand_regions/main'



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
    all_bedfile
    nonprot_bedfile
    prot_bedfile
    synonymous_bedfile
    exons_bedfile

    main:

    DNA2PROTEINMAPPING(mutations, exons_consensus_panel, all_samples_depths)


    // Create a channel for the domains file if omega_autodomains is true
    domains_ch = params.omega_autodomains ? domains_file : []  // .map{ it -> it[1]} : []
    exons_ch = params.omega_autoexons ? DNA2PROTEINMAPPING.out.panel_exons_bed.map{ it -> it[1]} : []

    // Create a channel for the hotspots bedfile if provided
    subgenic_ch = params.omega_subgenic_bedfile ? file(params.omega_subgenic_bedfile) : []

    if (params.omega_withingene){
        EXPANDREGIONSALL(all_consensus_panel, domains_ch, exons_ch, subgenic_ch)
        all_expanded_panel = EXPANDREGIONSALL.out.panel_increased.first()
        all_json_hotspots = EXPANDREGIONSALL.out.new_regions_json.first()

        EXPANDREGIONSNONPROT(nonprot_consensus_panel, domains_ch, exons_ch, subgenic_ch)
        nonprot_expanded_panel = EXPANDREGIONSNONPROT.out.panel_increased.first()
        nonprot_json_hotspots = EXPANDREGIONSNONPROT.out.new_regions_json.first()

        EXPANDREGIONSPROT(prot_consensus_panel, domains_ch, exons_ch, subgenic_ch)
        prot_expanded_panel = EXPANDREGIONSPROT.out.panel_increased.first()
        prot_json_hotspots = EXPANDREGIONSPROT.out.new_regions_json.first()

        EXPANDREGIONSSYNONYMOUS(synonymous_consensus_panel, domains_ch, exons_ch, subgenic_ch)
        synonymous_expanded_panel = EXPANDREGIONSSYNONYMOUS.out.panel_increased.first()
        synonymous_json_hotspots = EXPANDREGIONSSYNONYMOUS.out.new_regions_json.first()

        EXPANDREGIONSEXONS(exons_consensus_panel, domains_ch, exons_ch, subgenic_ch)
        exons_expanded_panel = EXPANDREGIONSEXONS.out.panel_increased.first()
        exons_json_hotspots = EXPANDREGIONSEXONS.out.new_regions_json.first()

    } else {
        all_expanded_panel = all_consensus_panel.first()
        nonprot_expanded_panel = nonprot_consensus_panel.first()
        prot_expanded_panel = prot_consensus_panel.first()
        synonymous_expanded_panel = synonymous_consensus_panel.first()
        exons_expanded_panel = exons_consensus_panel.first()

        all_json_hotspots = all_bedfile.first()
        nonprot_json_hotspots = nonprot_bedfile.first()
        prot_json_hotspots = prot_bedfile.first()
        synonymous_json_hotspots = synonymous_bedfile.first()
        exons_json_hotspots = exons_bedfile.first()
    }

    emit:

    all_consensus_expanded_panel        = all_expanded_panel
    nonprot_consensus_expanded_panel    = nonprot_expanded_panel
    prot_consensus_expanded_panel       = prot_expanded_panel
    synonymous_consensus_expanded_panel = synonymous_expanded_panel
    exons_consensus_expanded_panel      = exons_expanded_panel

    all_json_hotspots                   = all_json_hotspots
    nonprot_json_hotspots               = nonprot_json_hotspots
    prot_json_hotspots                  = prot_json_hotspots
    synonymous_json_hotspots            = synonymous_json_hotspots
    exons_json_hotspots                 = exons_json_hotspots

    dna2protein_mapping_depth_exons     = DNA2PROTEINMAPPING.out.depths_exons_positions.first()
    dna2protein_mapping_panel_exons     = DNA2PROTEINMAPPING.out.panel_exons_bed.first()
}
