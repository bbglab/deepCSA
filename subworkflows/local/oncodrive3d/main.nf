include { TABIX_BGZIPTABIX_QUERY    as SUBSETMUTATIONS      } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF                as SUBSETONCODRIVE3D   } from '../../../modules/local/subsetmaf/main'

include { ONCODRIVE3D_PREPROCESSING                         } from '../../../modules/local/bbgtools/oncodrive3d/preprocessing/main'
include { ONCODRIVE3D_RUN                                   } from '../../../modules/local/bbgtools/oncodrive3d/run/main'
include { ONCODRIVE3D_PLOT                                  } from '../../../modules/local/bbgtools/oncodrive3d/plot/main'
include { ONCODRIVE3D_PLOT_CHIMERAX                         } from '../../../modules/local/bbgtools/oncodrive3d/plot_chimerax/main'

// include { ONCODRIVE3D_COMP_PLOT } from '../../../modules/local/bbgtools/oncodrive3d/comparative_plot/main'



workflow ONCODRIVE3D_ANALYSIS{

    take:
    mutations
    mutabilities
    bedfile
    datasets
    annotations
    raw_vep

    main:

    SUBSETMUTATIONS(mutations, bedfile)
    SUBSETONCODRIVE3D(SUBSETMUTATIONS.out.subset)


    // mutations preprocessing
    if (params.o3d_raw_vep){
        ONCODRIVE3D_PREPROCESSING(SUBSETONCODRIVE3D.out.mutations, raw_vep)
        ONCODRIVE3D_PREPROCESSING.out.vep_output4o3d
        .join(mutabilities)
        .set{ muts_n_mutability }

    } else {
        SUBSETONCODRIVE3D.out.mutations
        .join(mutabilities)
        .set{ muts_n_mutability }

    }



    // if (datasets) {
    //      do something
    // } else {
    //     BUILDDATASETS
    // }


    ONCODRIVE3D_RUN(muts_n_mutability, datasets)

    if (params.o3d_plot) {
        ONCODRIVE3D_PLOT(ONCODRIVE3D_RUN.out.csv_genes
                            | join(ONCODRIVE3D_RUN.out.csv_pos)
                            | join(ONCODRIVE3D_RUN.out.mut_processed)
                            | join(ONCODRIVE3D_RUN.out.prob_processed)
                            | join(ONCODRIVE3D_RUN.out.seq_processed),
                            datasets,
                            annotations)
    }

    if (params.o3d_plot_chimerax) {
        ONCODRIVE3D_PLOT_CHIMERAX(ONCODRIVE3D_RUN.out.csv_genes
                                    | join(ONCODRIVE3D_RUN.out.csv_pos)
                                    | join(ONCODRIVE3D_RUN.out.prob_processed)
                                    | join(ONCODRIVE3D_RUN.out.seq_processed),
                                    datasets)
    }


    emit:

    results             = ONCODRIVE3D_RUN.out.csv_genes
    results_pos         = ONCODRIVE3D_RUN.out.csv_pos
    // files_mut           = ONCODRIVE3D_RUN.out.mut_processed
    // files_prob          = ONCODRIVE3D_RUN.out.prob_processed
    // files_seq           = ONCODRIVE3D_RUN.out.seq_processed
    // genes_plot          = ONCODRIVE3D_PLOT.out.genes_plot
    // summary_plot        = ONCODRIVE3D_PLOT.out.logodds_plot
    // volcano_plot        = ONCODRIVE3D_PLOT.out.volcano_plot
    // volcano_plot_gene   = ONCODRIVE3D_PLOT.out.volcano_plot_gene
    // results_pos_annot   = ONCODRIVE3D_PLOT.out.pos_annotated_csv
    // chimerax_defattr    = ONCODRIVE3D_PLOT_CHIMERAX.out.chimerax_defattr
    // chimerax_plot       = ONCODRIVE3D_PLOT_CHIMERAX.out.chimerax_plot
    // log_files           = ONCODRIVE3D_RUN.out.log.mix(ONCODRIVE3D_PLOT.out.log, ONCODRIVE3D_PLOT_CHIMERAX.out.log)

}
