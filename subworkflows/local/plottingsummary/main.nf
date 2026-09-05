


include { PLOT_SELECTION_METRICS            as PLOTSELECTION                    } from '../../../modules/local/plot/selection_metrics/main'
include { PLOT_SATURATION                   as PLOTSATURATION                   } from '../../../modules/local/plot/saturation/main'
include { PLOT_SATURATION_PROPORTIONS       as PLOTSATURATIONPROPORTIONS        } from '../../../modules/local/plot/saturation/proportions/main'
include { PLOT_INTERINDIVIDUAL_VARIABILITY  as PLOTINTERINDIVIDUALVARIABILITY   } from '../../../modules/local/plot/interindividual_variability/main'
include { COMPUTE_SATURATION_KINETICS       as COMPUTESATURATIONKINETICS        } from '../../../modules/local/saturation_kinetics/compute/main'
include { PLOT_SATURATION_KINETICS          as PLOTSTATURATIONKINETICS          } from '../../../modules/local/saturation_kinetics/plot/main'


workflow PLOTTING_SUMMARY {

    take:
    positive_selection_results_ready
    all_mutations
    all_mutdensities
    all_mutdensities_adjusted
    site_comparison
    all_samples_depth
    samples
    all_groups
    panel
    expanded_panel
    full_panel_rich
    seqinfo_df
    domain_df
    exons_depths_df
    groups_channel
    depths_indv
    relative_mutability
    omega_mutabilities


    main:

    pdb_tool_df   = params.annotations3d
                            ? channel.fromPath( "${params.annotations3d}/pdb_tool_df.tsv", checkIfExists: true).first()
                            : channel.empty()

    // think if we want to include this here
    // PLOTNEEDLES(muts_all_samples, sequence_information_df)


    if ( params.plot_only_allsamples ) {
        // plotting only for the entire cohort group
        channel.of([ [ id: "all_samples" ] ])
        .join( positive_selection_results_ready )
        .set{ groups_results }

        channel.of([ [ id: "all_samples" ] ])
        .join( all_mutations )
        .set{ groups_mutations }

    } else {
        // plotting for all groups
        positive_selection_results_ready
        .map { mut -> tuple(mut[0].id, mut) }
        .join(groups_channel)
        .map { it -> it[1] }
        .set { groups_results }

        all_mutations
        .map { mut -> tuple(mut[0].id, mut) }
        .join(groups_channel)
        .map { it -> it[1] }
        .set { groups_mutations }
    }

    groups_results
    .join( site_comparison )
    .set{ groups_results_sites }

    PLOTSELECTION(groups_results, seqinfo_df)
    // needles with consequence type
    // plot selection at cohort/group level, all the different methods available
    // plot selection per domain at cohort level

    PLOTSATURATION(groups_results_sites, all_samples_depth, panel, seqinfo_df, pdb_tool_df, domain_df, exons_depths_df)


    PLOTSATURATIONPROPORTIONS(groups_mutations, panel, full_panel_rich, expanded_panel)
    // plot gene + site selection
    // omega selection per domain in gene
    
    // plot saturation kinetics curves
    groups_mutations
    .join(depths_indv)
    .join(omega_mutabilities)
    .join(relative_mutability)
    .set{ groups_mutations_depths_n_mutability }
    COMPUTESATURATIONKINETICS(groups_mutations_depths_n_mutability, full_panel_rich, expanded_panel)
    // PLOTSTATURATIONKINETICS(groups_mutations, panel, full_panel_rich, expanded_panel)


    PLOTINTERINDIVIDUALVARIABILITY(samples, all_groups, panel,  all_mutdensities, all_mutdensities_adjusted)
    // heatmaps:
    //     mutations per gene/sample (total, SNV only, INDEL only, per consequence type)
    //     driver mutations per gene/sample
    //     mutation density
    //     mutation burden
    //     omega
    //     siteselection per group
    // define features based on PLOTMETRICSQC()


    emit:
    selection_plots     = PLOTSELECTION.out.plots
    saturation_plots    = PLOTSATURATION.out.plots

}
