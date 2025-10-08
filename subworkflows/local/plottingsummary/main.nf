


include { PLOT_SELECTION_METRICS            as PLOTSELECTION                    } from '../../../modules/local/plot/selection_metrics/main'
include { PLOT_SATURATION                   as PLOTSATURATION                   } from '../../../modules/local/plot/saturation/main'
include { PLOT_INTERINDIVIDUAL_VARIABILITY  as PLOTINTERINDIVIDUALVARIABILITY   } from '../../../modules/local/plot/interindividual_variability/main'



workflow PLOTTING_SUMMARY {

    take:
    positive_selection_results_ready
    all_mutdensities
    site_comparison
    all_samples_depth
    samples
    all_groups
    panel
    full_panel_rich
    seqinfo_df
    domain_df
    exons_depths_df
    groups_channel


    main:

    pdb_tool_df   = params.annotations3d
                            ? Channel.fromPath( "${params.annotations3d}/pdb_tool_df.tsv", checkIfExists: true).first()
                            : Channel.empty()

    // think if we want to include this here
    // PLOTNEEDLES(muts_all_samples, sequence_information_df)


    if ( params.plot_only_allsamples ) {
        // plotting only for the entire cohort group
        Channel.of([ [ id: "all_samples" ] ])
        .join( positive_selection_results_ready )
        .set{ groups_results }
    } else {
        // plotting for all groups
        positive_selection_results_ready
        .map { mut -> tuple(mut[0].id, mut) }
        .join(groups_channel)
        .map { it[1] }
        .set { groups_results }
    }

    groups_results
    .join( site_comparison )
    .set{ groups_results_sites }

    PLOTSELECTION(groups_results, seqinfo_df)
    // needles with consequence type
    // plot selection at cohort/group level, all the different methods available
    // plot selection per domain at cohort level

    PLOTSATURATION(groups_results_sites, all_samples_depth, panel, seqinfo_df, pdb_tool_df, domain_df, exons_depths_df)
    // plot gene + site selection
    // omega selection per domain in gene
    // ? plot saturation kinetics curves


    PLOTINTERINDIVIDUALVARIABILITY(samples, all_groups, panel,  all_mutdensities)
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
