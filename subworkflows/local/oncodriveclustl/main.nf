include { CREATECUSTOMBEDFILE    as ONCODRIVECLUSTLBED     } from '../../../modules/local/createpanels/custombedfile/main'

include { SUBSET_MAF             as SUBSET_ONCODRIVECLUSTL } from '../../../modules/local/subsetmaf/main'

include { ONCODRIVECLUSTL                                  } from '../../../modules/local/bbgtools/oncodriveclustl/main'
// include { ONCODRIVECLUSTL           as ONCODRIVECLUSTLSNVS    } from '../../../modules/local/bbgtools/oncodriveclustl/main'



workflow ONCODRIVECLUSTL_ANALYSIS{

    take:
    mutations
    mutabilities
    panel_file

    main:
    ch_versions = Channel.empty()

    ONCODRIVECLUSTLBED(panel_file)
    ch_versions = ch_versions.mix(ONCODRIVECLUSTLBED.out.versions)

    SUBSET_ONCODRIVECLUSTL(mutations)
    ch_versions = ch_versions.mix(ONCODRIVECLUSTLBED.out.versions)

    SUBSET_ONCODRIVECLUSTL.out.mutations
    .join(mutabilities)
    .set{ muts_n_mutability }

    ONCODRIVECLUSTL(muts_n_mutability,  ONCODRIVECLUSTLBED.out.bed)
    // ONCODRIVECLUSTLSNVS(muts_n_mutability,  ONCODRIVECLUSTLBED.out.bed)
    ch_versions = ch_versions.mix(ONCODRIVECLUSTL.out.versions)

    emit:
    results         = ONCODRIVECLUSTL.out.tsv          // channel: [ val(meta), file(results) ]
    // results_snvs    = ONCODRIVECLUSTLSNVS.out.tsv      // channel: [ val(meta), file(results) ]
    versions        = ch_versions                   // channel: [ versions.yml ]
}
