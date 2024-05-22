include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS          } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF    as SUBSET_ONCODRIVE3D } from '../../../modules/local/subsetmaf/main'

include { ONCODRIVE3D_RUN } from '../../../modules/local/bbgtools/oncodrive3d/run/main'

include { ONCODRIVE3D_PLOT } from '../../../modules/local/bbgtools/oncodrive3d/plot/main'

// include { ONCODRIVE3D_COMP_PLOT } from '../../../modules/local/bbgtools/oncodrive3d/comparative_plot/main'



workflow ONCODRIVE3D_ANALYSIS{

    take:
    mutations
    mutabilities
    bedfile
    datasets
    annotations

    main:
    ch_versions = Channel.empty()

    SUBSETMUTATIONS(mutations, bedfile)
    ch_versions = ch_versions.mix(SUBSETMUTATIONS.out.versions)

    SUBSET_ONCODRIVE3D(SUBSETMUTATIONS.out.subset)
    ch_versions = ch_versions.mix(SUBSET_ONCODRIVE3D.out.versions)

    SUBSET_ONCODRIVE3D.out.mutations
    .join(mutabilities)
    .set{ muts_n_mutability }

    // if (datasets) {
    //      do something
    // } else {
    //     BUILDDATASETS
    // }

    ONCODRIVE3D_RUN(muts_n_mutability, datasets)
    ch_versions = ch_versions.mix(ONCODRIVE3D_RUN.out.versions)

    ONCODRIVE3D_PLOT(ONCODRIVE3D_RUN.out.csv_genes,
                     ONCODRIVE3D_RUN.out.csv_pos,
                     ONCODRIVE3D_RUN.out.mut_processed,
                     ONCODRIVE3D_RUN.out.prob_processed,
                     ONCODRIVE3D_RUN.out.seq_processed,
                     datasets,
                     annotations)
    ch_versions = ch_versions.mix(ONCODRIVE3D_PLOT.out.versions)

    emit:
    results           = ONCODRIVE3D_RUN.out.csv_genes
    results_pos       = ONCODRIVE3D_RUN.out.csv_pos
    files_mut         = ONCODRIVE3D_RUN.out.mut_processed
    files_prob        = ONCODRIVE3D_RUN.out.prob_processed
    files_seq         = ONCODRIVE3D_RUN.out.seq_processed
    genes_plot        = ONCODRIVE3D_PLOT.out.genes_plot
    summary_plot      = ONCODRIVE3D_PLOT.out.summary_plot
    log_files         = ONCODRIVE3D_RUN.out.log.mix(ONCODRIVE3D_PLOT.out.log)
    

    versions = ch_versions           // channel: [ versions.yml ]

}
