
include { PLOT_MUTDENSITY_QC    as PLOTMUTDENSITYQC     } from '../../../modules/local/plot/qc/mutation_densities/main'
include { PLOT_METRICS_VS_DEPTH_QC as PLOTMETRICSVSDEPTHQC } from '../../../modules/local/plot/qc/metrics_vs_depth/main'
include { ANNOTATE_OMEGA_QC     as APPLYOMEGAQC         } from '../../../modules/local/plot/qc/annotate_omega/main'
include { PLOT_MUTATION_SPECIFIC       as PLOTMUTATIONSPECIFIC  } from '../../../modules/local/plot/qc/mutation_specific/main'


workflow PLOTTING_QC {

    take:
    all_mutations
    // positive_selection_results_ready
    all_mutdensities
    all_adjusted_mutdensities
    all_omegas_globalloc
    average_depth_gene_sample
    // all_samples_depth
    // all_groups
    all_omegas
    panel
    groups_definition
    group_name
    // full_panel_rich
    // seqinfo_df
    // domain_df
    // exons_depths_df


    main:

    // Channel.of([ [ id: "all_samples" ] ])
    // .join( all_mutations )
    // .set{ mutations }
    PLOTMUTATIONSPECIFIC(all_mutations)

    // plotting only for the entire cohort group
    // channel.of([ [ id: "all_samples" ] ])
    // .join( positive_selection_results_ready )
    // .set{ all_samples_results }

    PLOTMUTDENSITYQC(all_mutdensities, panel, groups_definition, group_name)

    PLOTMETRICSVSDEPTHQC(
        all_mutdensities,
        // average_depth_gene_sample channel emits tuples [meta, depth_per_gene_per_sample.tsv];
        // pass only the file path to the plotting module.
        average_depth_gene_sample.map { it -> it[1] },
        groups_definition,
        'all_samples', // group_name, // only activate all samples while developing; eventually will want to loop over all groups
        all_adjusted_mutdensities,
        all_omegas_globalloc
    )
    // mutation density per gene cohort-level
    // mutation density per gene & sample
    //      synonymous
    //      non-protein-affecting
    // pending:
    //      protein-affecting
    //          truncating
    //          missense

    APPLYOMEGAQC(all_omegas, PLOTMUTDENSITYQC.out.compiled_flagged.collect())

    emit:
    mutdensity_plots        = PLOTMUTDENSITYQC.out.plots
    metrics_vs_depth_plots  = PLOTMETRICSVSDEPTHQC.out.plots
    metrics_vs_depth_tables = PLOTMETRICSVSDEPTHQC.out.tables
    flagged_omegas          = APPLYOMEGAQC.out.all_omegas_annotated

}
