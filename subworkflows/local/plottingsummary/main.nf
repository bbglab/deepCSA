


include { PLOT_SELECTION_METRICS            as PLOTSELECTION                    } from '../../../modules/local/plot/selection_metrics/main'
include { PLOT_SATURATION                   as PLOTSATURATION                   } from '../../../modules/local/plot/saturation/main'
// include { PLOT_INTERINDIVIDUAL_VARIABILITY  as PLOTINTERINDIVIDUALVARIABILITY   } from '../../../modules/local/plot/selection_metrics/main'
// include { PLOT_METRICS_QC                   as PLOTMETRICSQC                    } from '../../../modules/local/plot/selection_metrics/main'



workflow PLOTTING_SUMMARY {

    take:
    positive_selection_results_ready
    all_mutrates
    samples
    all_groups
    panel
    full_panel_rich
    seqinfo_df
    all_samples_depth


    main:

    pdb_tool_df   = params.annotations3d
                            ? Channel.fromPath( "${params.annotations3d}/pdb_tool_df.tsv", checkIfExists: true).first()
                            : Channel.empty()

    // replace by domains file
    domain_df      = params.annotations3d
                            ? Channel.fromPath( "${params.annotations3d}/uniprot_feat.tsv", checkIfExists: true).first()
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



    PLOTSATURATION(all_samples_results, all_samples_depth, panel, seqinfo_df, pdb_tool_df, domain_df)
    // // plot gene + site selection + omega selection per domain in gene
    // // ? plot saturation kinetics curves



    // PLOTMETRICSQC(all_mutrates, )
    // // mutation density per gene cohort-level
    // // mutation density per gene & sample
    // //      synonymous
    // //      protein-affecting
    // //          truncating
    // //          missense


    // PLOTINTERINDIVIDUALVARIABILITY()
    // // heatmap driver mutations per gene/sample
    // // other heatmaps:
    // //     mutation density
    // //     mutation burden
    // //     omega
    // //     siteselection per group
    // // define features based on PLOTMETRICSQC()


    emit:
    selection_plots = PLOTSELECTION.out.plots

}
