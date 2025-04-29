include { CREATECUSTOMBEDFILE    as ONCODRIVECLUSTLBED     } from '../../../modules/local/createpanels/custombedfile/main'

include { SUBSET_MAF             as SUBSETONCODRIVECLUSTL } from '../../../modules/local/subsetmaf/main'

include { ONCODRIVECLUSTL                                  } from '../../../modules/local/bbgtools/oncodriveclustl/main'
// include { ONCODRIVECLUSTL           as ONCODRIVECLUSTLSNVS    } from '../../../modules/local/bbgtools/oncodriveclustl/main'



workflow ONCODRIVECLUSTL_ANALYSIS{

    take:
    mutations
    mutabilities
    panel_file

    main:

    ONCODRIVECLUSTLBED(panel_file)

    SUBSETONCODRIVECLUSTL(mutations)

    SUBSETONCODRIVECLUSTL.out.mutations
    .join(mutabilities)
    .set{ muts_n_mutability }

    ONCODRIVECLUSTL(muts_n_mutability,  ONCODRIVECLUSTLBED.out.bed)
    // ONCODRIVECLUSTLSNVS(muts_n_mutability,  ONCODRIVECLUSTLBED.out.bed)

    emit:
    results         = ONCODRIVECLUSTL.out.tsv          // channel: [ val(meta), file(results) ]
    // results_snvs    = ONCODRIVECLUSTLSNVS.out.tsv      // channel: [ val(meta), file(results) ]
}
