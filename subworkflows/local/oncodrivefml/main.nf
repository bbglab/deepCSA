include { CREATECUSTOMBEDFILE    as ONCODRIVEFMLBED     } from '../../../modules/local/createpanels/custombedfile/main'

include { SUBSET_MAF             as SUBSET_ONCODRIVEFML } from '../../../modules/local/subsetmaf/main'

include { ONCODRIVEFML                                  } from '../../../modules/local/bbgtools/oncodrivefml/main'
include { ONCODRIVEFML           as ONCODRIVEFMLSNVS    } from '../../../modules/local/bbgtools/oncodrivefml/main'



workflow ONCODRIVEFML_ANALYSIS{

    take:
    mutations
    mutabilities
    panel_file
    cadd_scores

    main:
    ch_versions = Channel.empty()

    ONCODRIVEFMLBED(panel_file)
    ch_versions = ch_versions.mix(ONCODRIVEFMLBED.out.versions)

    SUBSET_ONCODRIVEFML(mutations)
    ch_versions = ch_versions.mix(SUBSET_ONCODRIVEFML.out.versions)

    SUBSET_ONCODRIVEFML.out.mutations
    .join( mutabilities )
    .set{ muts_n_mutability }

    ONCODRIVEFML(muts_n_mutability,  ONCODRIVEFMLBED.out.bed, cadd_scores)
    ONCODRIVEFMLSNVS(muts_n_mutability,  ONCODRIVEFMLBED.out.bed, cadd_scores)
    ch_versions = ch_versions.mix(ONCODRIVEFML.out.versions)

    emit:
    results         = ONCODRIVEFML.out.tsv          // channel: [ val(meta), file(results) ]
    results_snvs    = ONCODRIVEFMLSNVS.out.tsv      // channel: [ val(meta), file(results) ]
    versions        = ch_versions                   // channel: [ versions.yml ]
}
