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
    ch_versions = Channel.empty()

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETDEPTHS(depth, bedfile)
    ch_versions = ch_versions.mix(SUBSETDEPTHS.out.versions)

    ONCODRIVEFMLBED(panel)
    ch_versions = ch_versions.mix(ONCODRIVEFMLBED.out.versions)

    PLOTDEPTHS(SUBSETDEPTHS.out.subset, ONCODRIVEFMLBED.out.bed)
    ch_versions = ch_versions.mix(PLOTDEPTHS.out.versions)


    emit:
    plots = PLOTDEPTHS.out.plots
    versions = ch_versions                // channel: [ versions.yml ]
}
