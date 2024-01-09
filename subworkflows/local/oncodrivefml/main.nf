include { SUBSET_MAF    as SUBSET_ONCODRIVEFML } from '../../../modules/local/subsetmaf/main'
include { ONCODRIVEFML } from '../../../modules/local/bbgtools/oncodrivefml/main'



workflow ONCODRIVEFML_ANALYSIS{

    take:
    // inputs
    muts
    bedfile

    main:
    ch_versions = Channel.empty()


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

    SUBSET_ONCODRIVEFML.out.mutations
    .join(meta_mutabs)
    .set{ muts_n_mutability }

    ONCODRIVEFML(muts_n_mutability,    bedfile)
    ch_versions = ch_versions.mix(ONCODRIVEFML.out.versions)

    // TODO: add a second run of ONCODRIVEFML only with SNVs --no-indels

    emit:
    results = ONCODRIVEFML.out.tsv   // channel: [ val(meta), file(results) ]
    versions = ch_versions           // channel: [ versions.yml ]
}
