include { TABIX_BGZIPTABIX_QUERY    as QUERYMUTATIONS           } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { SUBSET_MAF                as SUBSETMUTDENSITYADJUSTED } from '../../../modules/local/subsetmaf/main'

include { MUTATION_DENSITY          as MUTDENSITYADJ            } from '../../../modules/local/mut_density/adjusted/main'


workflow MUTATION_DENSITY {
    take:
    mutations
    depth
    bedfile
    panel
    mutability_file
    trinucleotide_counts

    main:

    // Intersect BED of all sites of interest with mutations
    QUERYMUTATIONS(mutations, bedfile)

    SUBSETMUTDENSITYADJUSTED(QUERYMUTATIONS.out.subset)

    SUBSETMUTDENSITYADJUSTED.out.mutations
    .join(depth)
    .join(mutability_file)
    .set{ mutations_n_depth_n_profile }

    MUTDENSITYADJ(mutations_n_depth_n_profile, panel, trinucleotide_counts)

    // Concatenate all outputs into a single file
    MUTDENSITYADJ.out.mutdensities.map{ it -> it[1]}.flatten()
    .set{ all_adjusted_mutdensities }
    all_adjusted_mutdensities.collectFile(name: "all_adjusted_mutdensities.tsv", storeDir:"${params.outdir}/mutdensity_adjusted", skip: 1, keepHeader: true).first().set{ all_adjusted_mutdensities_file }

    MUTDENSITYADJ.out.mutdensities_flat.map{ it -> it[1]}.flatten()
    .set{ all_adjusted_mutdensities_flat }
    all_adjusted_mutdensities_flat.collectFile(name: "all_adjusted_mutdensities_flat.tsv", storeDir:"${params.outdir}/mutdensity_adjusted", skip: 1, keepHeader: true)

    channel.of([ [ id: "all_samples" ] ])
    .join( MUTDENSITYADJ.out.mutdensities )
    .set{ all_samples_adj_mutdensity }



    emit:
    mutdensities            = MUTDENSITYADJ.out.mutdensities
    mutdensities_flat       = MUTDENSITYADJ.out.mutdensities_flat
    mutdensities_plots      = MUTDENSITYADJ.out.mutdensities_plots
    all_mutdensities        = all_adjusted_mutdensities
    all_mutdensities_flat   = all_adjusted_mutdensities_flat
    all_samples_mutdensity  = all_samples_adj_mutdensity
    all_adj_mutdens_file    = all_adjusted_mutdensities_file
}
