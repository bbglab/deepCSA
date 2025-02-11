include { CREATECUSTOMBEDFILE    as ONCODRIVEFMLBED     } from '../../../modules/local/createpanels/custombedfile/main'

include { SUBSET_MAF             as SUBSETONCODRIVEFML  } from '../../../modules/local/subsetmaf/main'

include { ONCODRIVEFML                                  } from '../../../modules/local/bbgtools/oncodrivefml/main'
include { ONCODRIVEFML           as ONCODRIVEFMLSNVS    } from '../../../modules/local/bbgtools/oncodrivefml/main'



workflow ONCODRIVEFML_ANALYSIS{

    take:
    mutations
    mutabilities
    panel_file
    cadd_scores
    mode

    main:
    ONCODRIVEFMLBED(panel_file)

    SUBSETONCODRIVEFML(mutations)

    SUBSETONCODRIVEFML.out.mutations
    .join( mutabilities )
    .set{ muts_n_mutability }


    ONCODRIVEFML(muts_n_mutability,  ONCODRIVEFMLBED.out.bed, cadd_scores, mode)
    ONCODRIVEFMLSNVS(muts_n_mutability,  ONCODRIVEFMLBED.out.bed, cadd_scores, mode)


    emit:
    results                = ONCODRIVEFML.out.tsv             // channel: [ val(meta), file(results) ]
    results_snvs           = ONCODRIVEFMLSNVS.out.tsv         // channel: [ val(meta), file(results) ]
    results_folder         = ONCODRIVEFML.out.folder          // channel: [ val(meta), file(results) ]
    results_snvs_folder    = ONCODRIVEFMLSNVS.out.folder      // channel: [ val(meta), file(results) ]

}
