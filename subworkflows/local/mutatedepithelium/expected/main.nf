
include { TABIX_BGZIPTABIX_QUERY        as SUBSETDEPTHS      } from '../../../../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY        as SUBSETMUTATIONS   } from '../../../../modules/nf-core/tabix/bgziptabixquery/main'
include { SUBSET_MAF                    as SUBSETMUTEPIVAF  } from '../../../../modules/local/subsetmaf/main'
include { EXP_MUTRATE                   as EXPMUTRATE        } from '../../../../modules/local/expected_mutrate/main'
include { CREATECUSTOMBEDFILE           as INTERVALSBED      } from '../../../../modules/local/createpanels/custombedfile/main'

workflow EXPECTED_MUTRATE {

    take:
    mutations
    bedfile
    panel
    depth
    raw_annotation

    main:

    SUBSETMUTATIONS(mutations, bedfile)

    // Intersect BED of all sites with BED of sample filtered sites
    SUBSETDEPTHS(depth, bedfile)

    // SUBSET_MUTEPIVAF(SUBSETMUTATIONS.out.subset)
    INTERVALSBED(panel)

    features_table = params.features_table ? Channel.fromPath( params.features_table, checkIfExists: true) : Channel.fromPath(params.input)
    EXPMUTRATE(panel,
                SUBSETMUTATIONS.out.subset,
                SUBSETDEPTHS.out.subset,
                raw_annotation,
                INTERVALSBED.out.bed,
                features_table
                )

    emit:
    // TODO add some other output
    // mut_epi_sample  = COMPUTEMUTATEDEPITHELIUM.out.mutated_epi_sample
    refcds_object       = EXPMUTRATE.out.rds_file.first()
    refcds_object_rda   = EXPMUTRATE.out.rda_file.first()

}
