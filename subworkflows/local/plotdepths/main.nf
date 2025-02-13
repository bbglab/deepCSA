include { TABIX_BGZIPTABIX_QUERY    as SUBSETDEPTHS        } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { PLOT_DEPTHS               as PLOTDEPTHS          } from '../../../modules/local/plot/depths_summary/main'

include { MUTRATE                   as MUTRATE             } from '../../../modules/local/computemutrate/main'

include { CREATECUSTOMBEDFILE       as ONCODRIVEFMLBED     } from '../../../modules/local/createpanels/custombedfile/main'

workflow PLOT_DEPTHS {

    take:
    depth
    bedfile
    panel

    main:

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETDEPTHS(depth, bedfile)

    ONCODRIVEFMLBED(panel)

    PLOTDEPTHS(SUBSETDEPTHS.out.subset, ONCODRIVEFMLBED.out.bed)


    emit:
    plots = PLOTDEPTHS.out.plots
}
