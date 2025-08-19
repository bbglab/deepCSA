


include { PLOT_SELECTION_METRICS            as PLOTSELECTION                    } from '../../../modules/local/plot/selection_metrics/main'
include { PLOT_SATURATION                   as PLOTSATURATION                   } from '../../../modules/local/plot/saturation/main'
// include { PLOT_INTERINDIVIDUAL_VARIABILITY  as PLOTINTERINDIVIDUALVARIABILITY   } from '../../../modules/local/plot/selection_metrics/main'
// include { PLOT_METRICS_QC                   as PLOTMETRICSQC                    } from '../../../modules/local/plot/selection_metrics/main'



workflow PLOTTING_SUMMARY {

    take:
    positive_selection_results_ready
    all_mutrates
    site_comparison
    all_samples_depth
    samples
    all_groups
    panel
    full_panel_rich
    seqinfo_df
    domain_df
    exons_depths_df


    main:

    pdb_tool_df   = params.annotations3d
                            ? Channel.fromPath( "${params.annotations3d}/pdb_tool_df.tsv", checkIfExists: true).first()
                            : Channel.empty()

    // think if we want to include this here
    // PLOTNEEDLES(muts_all_samples, sequence_information_df)


    // plotting only for the entire cohort group
    Channel.of([ [ id: "all_samples" ] ])
    .join( positive_selection_results_ready )
    .set{ all_samples_results }

    PLOTSELECTION(all_samples_results, seqinfo_df)
    // needles with consequence type
    // plot selection at cohort/group level, all the different methods available
    // plot selection per domain at cohort level

    Channel.of([ [ id: "all_samples" ] ])
    .join( site_comparison )
    .set{ all_samples_sites }

    PLOTSATURATION(all_samples_results, all_samples_depth, all_samples_sites, panel, seqinfo_df, pdb_tool_df, domain_df, exons_depths_df)
    // plot gene + site selection
    // omega selection per domain in gene
    // ? plot saturation kinetics curves



    // PLOTMETRICSQC(all_mutrates, )
    // // mutation density per gene cohort-level
    // // mutation density per gene & sample
    // //      synonymous
    // //      protein-affecting
    // //          truncating
    // //          missense


    // PLOTINTERINDIVIDUALVARIABILITY(all_mutrates, all_samples_depth, )
    // // heatmaps:
    // //     mutations per gene/sample (total, SNV only, INDEL only, per consequence type)
    // //     driver mutations per gene/sample
    // //     mutation density
    // //     mutation burden
    // //     omega
    // //     siteselection per group
    // // define features based on PLOTMETRICSQC()


    emit:
    selection_plots     = PLOTSELECTION.out.plots
    saturation_plots    = PLOTSATURATION.out.plots

}
