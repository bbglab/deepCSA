include { CREATECUSTOMBEDFILE    as ONCODRIVEFMLBED     } from '../../../modules/local/createpanels/custombedfile/main'

include { SUBSET_MAF             as SUBSET_ONCODRIVEFML } from '../../../modules/local/subsetmaf/main'

include { ONCODRIVEFML                                  } from '../../../modules/local/bbgtools/oncodrivefml/main'
include { ONCODRIVEFML           as ONCODRIVEFMLSNVS    } from '../../../modules/local/bbgtools/oncodrivefml/main'



workflow ONCODRIVEFML_ANALYSIS{

    take:
    // inputs
    muts
    panel_file

    main:
    ch_versions = Channel.empty()


    ONCODRIVEFMLBED(panel_file)
    ch_versions = ch_versions.mix(ONCODRIVEFMLBED.out.versions)

    //
    // Separate input mutations and mutabilities
    //
    muts.
    map{ it -> [it[0], it[1]]}.
    set{ meta_muts }

    muts.
    map{ it -> [it[0], it[2], it[3]]}.
    set{ meta_mutabs }

    SUBSET_ONCODRIVEFML(meta_muts)
    ch_versions = ch_versions.mix(SUBSET_ONCODRIVEFML.out.versions)

    SUBSET_ONCODRIVEFML.out.mutations
    .join(meta_mutabs)
    .set{ muts_n_mutability }

    ONCODRIVEFML(muts_n_mutability,  ONCODRIVEFMLBED.out.bed)
    ONCODRIVEFMLSNVS(muts_n_mutability,  ONCODRIVEFMLBED.out.bed)
    ch_versions = ch_versions.mix(ONCODRIVEFML.out.versions)

    emit:
    results         = ONCODRIVEFML.out.tsv          // channel: [ val(meta), file(results) ]
    results_snvs    = ONCODRIVEFMLSNVS.out.tsv      // channel: [ val(meta), file(results) ]
    versions        = ch_versions                   // channel: [ versions.yml ]
}
