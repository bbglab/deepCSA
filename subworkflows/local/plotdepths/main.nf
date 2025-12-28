include { TABIX_BGZIPTABIX_QUERY    as QUERYDEPTHS     } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { PLOT_DEPTHS               as DEPTHSSUMMARY    } from '../../../modules/local/plot/depths_summary/main'

include { CREATECUSTOMBEDFILE       as CUSTOMBEDFILE    } from '../../../modules/local/createpanels/custombedfile/main'


workflow PLOT_DEPTHS {

    take:
    depth
    bedfile
    panel

    main:

    // Intersect BED of all sites with BED of sample filtered sites
    QUERYDEPTHS(depth, bedfile)

    CUSTOMBEDFILE(panel)

    DEPTHSSUMMARY(QUERYDEPTHS.out.subset, CUSTOMBEDFILE.out.bed)    


    emit:
    plots           = DEPTHSSUMMARY.out.plots
    average_depth   = DEPTHSSUMMARY.out.average_per_sample

}
