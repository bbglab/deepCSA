include { TABIX_BGZIPTABIX_QUERY_INDEX  as SUBSETPILEUP                 } from '../../../modules/nf-core/tabix/bgziptabixqueryindex/main'

include { COMPUTE_FRAGMENT_COORDS       as COMPUTEFRAGMENTSCOORDS       } from '../../../modules/local/fragments_from_bam/main'

include { CREATECUSTOMBEDFILE           as READSPOSBED                  } from '../../../modules/local/createpanels/custombedfile/main'
include { READS_PER_REGION              as READSPERREGION               } from '../../../modules/local/reads_per_region/main'

include { TABIX_BGZIPTABIX_QUERY        as SUBSETMUTATIONS              } from '../../../modules/nf-core/tabix/bgziptabixquery/main'
include { SUBSET_MAF                    as SUBSET_MUTEPI                } from '../../../modules/local/subsetmaf/main'

include { COMPUTE_MUTATED_EPITHELIUM    as COMPUTEMUTATEDEPITHELIUM     } from '../../../modules/local/computemutatedepithelium/main'




workflow MUTATED_EPITHELIUM {

    take:
    mutations
    bedfile
    panel
    bamfile
    pileup


    main:
    ch_versions = Channel.empty()

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETPILEUP(pileup, bedfile)
    ch_versions = ch_versions.mix(SUBSETPILEUP.out.versions)

    COMPUTEFRAGMENTSCOORDS(bamfile)
    ch_versions = ch_versions.mix(COMPUTEFRAGMENTSCOORDS.out.versions)

    READSPOSBED(panel)
    ch_versions = ch_versions.mix(READSPOSBED.out.versions)

    SUBSETPILEUP.out.subset
    .join(COMPUTEFRAGMENTSCOORDS.out.fragments)
    .set{ pileup_n_fragments }

    READSPERREGION(pileup_n_fragments, READSPOSBED.out.bed)
    ch_versions = ch_versions.mix(READSPERREGION.out.versions)

    // TODO
    // see if we should be using the bedfile we generated above here
    // or it does not matter
    SUBSETMUTATIONS(mutations, bedfile)
    ch_versions = ch_versions.mix(SUBSETMUTATIONS.out.versions)

    SUBSET_MUTEPI(SUBSETMUTATIONS.out.subset)
    ch_versions = ch_versions.mix(SUBSET_MUTEPI.out.versions)

    SUBSET_MUTEPI.out.mutations
    .join(READSPERREGION.out.read_counts.map{it -> [ ["id" : it[0].id] , it[1] ]  })
    .set{ mutations_n_reads }

    COMPUTEMUTATEDEPITHELIUM(mutations_n_reads)



    emit:
    read_counts     = READSPERREGION.out.read_counts
    mut_epi_exon    = COMPUTEMUTATEDEPITHELIUM.out.mutated_epi_exon
    mut_epi_gene    = COMPUTEMUTATEDEPITHELIUM.out.mutated_epi_gene
    mut_epi_sample  = COMPUTEMUTATEDEPITHELIUM.out.mutated_epi_sample
    versions        = ch_versions                // channel: [ versions.yml ]

}
