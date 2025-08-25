include { TABIX_BGZIPTABIX_QUERY as SUBSETDEPTHS } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { PLOT_DEPTHS as DEPTHSSUMMARY           } from '../../../modules/local/plot/depths_summary/main'

include { CREATECUSTOMBEDFILE as ONCODRIVEFMLBED } from '../../../modules/local/createpanels/custombedfile/main'


workflow PLOT_DEPTHS {
    take:
    depth
    bedfile
    panel

    main:

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETDEPTHS(depth, bedfile)

    ONCODRIVEFMLBED(panel)

    DEPTHSSUMMARY(SUBSETDEPTHS.out.subset, ONCODRIVEFMLBED.out.bed)

    emit:
    plots         = DEPTHSSUMMARY.out.plots
    average_depth = DEPTHSSUMMARY.out.average_per_sample
}
