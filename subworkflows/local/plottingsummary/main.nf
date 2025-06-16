


include { PLOT_SELECTION_METRICS            as PLOTSELECTION                    } from '../../../modules/local/plot/selection_metrics/main'
include { PLOT_SATURATION                   as PLOTSATURATION                   } from '../../../modules/local/plot/selection_metrics/main'
include { PLOT_INTERINDIVIDUAL_VARIABILITY  as PLOTINTERINDIVIDUALVARIABILITY   } from '../../../modules/local/plot/selection_metrics/main'
include { PLOT_METRICS_QC                   as PLOTMETRICSQC                    } from '../../../modules/local/plot/selection_metrics/main'



workflow PLOTTING_SUMMARY {

    take:
    positive_selection_results_ready
    all_mutrates
    samples
    groups
    all_groups
    panel
    full_panel_rich
    seqinfo_df

    main:

    // think if we want to include this here
    // PLOTNEEDLES(muts_all_samples, sequence_information_df)



    PLOTSELECTION(positive_selection_results_ready, seqinfo_df)
    // fig 2a omega per sample
    // needles with consequence type

    // plot selection at cohort/group level, all the different methods available

    // plot selection per domain at cohort level



    PLOTSATURATION()
    // plot gene + site selection + omega selection per domain in gene
    // ? plot saturation kinetics curves



    PLOTMETRICSQC(all_mutrates, )
    // mutation density per gene cohort-level
    // mutation density per gene & sample
    //      synonymous
    //      protein-affecting
    //          truncating
    //          missense


    PLOTINTERINDIVIDUALVARIABILITY()
    // heatmap driver mutations per gene/sample
    // other heatmaps:
    //     mutation density
    //     mutation burden
    //     omega
    //     siteselection per group
    // define features based on PLOTMETRICSQC()


    emit:
    selection_plots = PLOTSELECTION.out.plots

}
