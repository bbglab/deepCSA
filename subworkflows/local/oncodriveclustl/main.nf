include { CREATECUSTOMBEDFILE    as ONCODRIVECLUSTLBED     } from '../../../modules/local/createpanels/custombedfile/main'

include { SUBSET_MAF             as SUBSET_ONCODRIVECLUSTL } from '../../../modules/local/subsetmaf/main'

include { ONCODRIVECLUSTL                                  } from '../../../modules/local/bbgtools/oncodriveclustl/main'
// include { ONCODRIVECLUSTL           as ONCODRIVEFMLSNVS    } from '../../../modules/local/bbgtools/oncodrivefml/main'



workflow ONCODRIVECLUSTL_ANALYSIS{

    take:
    muts
    panel_file

    main:
    ch_versions = Channel.empty()


    ONCODRIVECLUSTLBED(panel_file)
    ch_versions = ch_versions.mix(ONCODRIVECLUSTLBED.out.versions)

    //
    // Separate input mutations and mutabilities
    //
    muts.
    map{ it -> [it[0], it[1]]}.
    set{ meta_muts }

    muts.
    map{ it -> [it[0], it[2], it[3]]}.
    set{ meta_mutabs }

    SUBSET_ONCODRIVECLUSTL(meta_muts)
    ch_versions = ch_versions.mix(ONCODRIVECLUSTLBED.out.versions)

    SUBSET_ONCODRIVECLUSTL.out.mutations
    .join(meta_mutabs)
    .set{ muts_n_mutability }

    ONCODRIVECLUSTL(muts_n_mutability,  ONCODRIVECLUSTLBED.out.bed)
    // ONCODRIVECLUSTLSNVS(muts_n_mutability,  ONCODRIVECLUSTLBED.out.bed)
    ch_versions = ch_versions.mix(ONCODRIVECLUSTL.out.versions)

    emit:
    results         = ONCODRIVECLUSTL.out.tsv          // channel: [ val(meta), file(results) ]
    // results_snvs    = ONCODRIVEFMLSNVS.out.tsv      // channel: [ val(meta), file(results) ]
    versions        = ch_versions                   // channel: [ versions.yml ]
}
