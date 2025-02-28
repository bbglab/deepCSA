
include { SAMTOOLS_MPILEUP              as PILEUPBAMALL                 } from '../../../../modules/nf-core/samtools/mpileup/main'
include { TABIX_BGZIPTABIX_QUERY_INDEX  as SUBSETPILEUP                 } from '../../../../modules/nf-core/tabix/bgziptabixqueryindex/main'

include { COMPUTE_FRAGMENT_COORDS       as COMPUTEFRAGMENTSCOORDS       } from '../../../../modules/local/mutated_cells_from_reads/fragments_from_bam/main'
include { READS_PER_REGION              as READSPERREGION               } from '../../../../modules/local/mutated_cells_from_reads/reads_per_region/main'

include { CREATECUSTOMBEDFILE           as READSPOSBED                  } from '../../../../modules/local/createpanels/custombedfile/main'

include { TABIX_BGZIPTABIX_QUERY        as SUBSETMUTATIONS              } from '../../../../modules/nf-core/tabix/bgziptabixquery/main'
include { SUBSET_MAF                    as SUBSETMUTEPI                 } from '../../../../modules/local/subsetmaf/main'
include { COMPUTE_MUTATED_EPITHELIUM    as COMPUTEMUTATEDEPITHELIUM     } from '../../../../modules/local/mutated_cells_from_reads/compute/main'




workflow MUTATED_CELLS_READS {

    take:
    mutations
    bedfile
    panel
    pileup_bam_bai
    fasta


    main:

    pileup_bam_bai.map{ it -> [ it[0], it[1] ]}
    .set{ bamfile }

    ch_bamall_bai_bed = pileup_bam_bai.combine(bedfile.map{ it -> [ it[1] ]}  )
    PILEUPBAMALL(ch_bamall_bai_bed, fasta)

    PILEUPBAMALL.out.mpileup.map{ it -> [it[0], it[1]] }
    .set{ pileup }

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETPILEUP(pileup, bedfile)

    COMPUTEFRAGMENTSCOORDS(bamfile)

    READSPOSBED(panel)

    SUBSETPILEUP.out.subset
    .join(COMPUTEFRAGMENTSCOORDS.out.fragments)
    .set{ pileup_n_fragments }

    READSPERREGION(pileup_n_fragments, READSPOSBED.out.bed)

    // TODO
    // see if we should be using the bedfile we generated above here
    // or it does not matter
    SUBSETMUTATIONS(mutations, bedfile)

    SUBSETMUTEPI(SUBSETMUTATIONS.out.subset)

    SUBSETMUTEPI.out.mutations
    .join(READSPERREGION.out.read_counts.map{it -> [ ["id" : it[0].id] , it[1] ]  })
    .set{ mutations_n_reads }

    COMPUTEMUTATEDEPITHELIUM(mutations_n_reads)


    emit:
    read_counts     = READSPERREGION.out.read_counts
    mut_epi_exon    = COMPUTEMUTATEDEPITHELIUM.out.mutated_epi_exon
    mut_epi_gene    = COMPUTEMUTATEDEPITHELIUM.out.mutated_epi_gene
    mut_epi_sample  = COMPUTEMUTATEDEPITHELIUM.out.mutated_epi_sample

}
