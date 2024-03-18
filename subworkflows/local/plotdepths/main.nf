include { TABIX_BGZIPTABIX_QUERY    as SUBSETDEPTHS             } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { PLOTDEPTHS                as PLOTDEPTHS               } from '../../../modules/local/plot/depths_summary/main'

include { MUTRATE as MUTRATE } from '../../../modules/local/computemutrate/main'


workflow PLOT_DEPTHS {

    take:
    depth
    bedfile
    panel

    main:
    ch_versions = Channel.empty()

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETDEPTHS(depth, bedfile)
    ch_versions = ch_versions.mix(SUBSETDEPTHS.out.versions)

    PLOTDEPTHS(SUBSETDEPTHS.out.subset, panel)
    ch_versions = ch_versions.mix(PLOTDEPTHS.out.versions)

    emit:
    plots = PLOTDEPTHS.out.plots
    versions = ch_versions                // channel: [ versions.yml ]
}
