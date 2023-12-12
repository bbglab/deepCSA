include { ONCODRIVEFML } from '../../../modules/local/bbgtools/oncodrivefml/main'

workflow ONCODRIVEFML_ANALYSIS{

    take:
    // inputs
    muts
    mutabs
    mutabind
    bedfile

    main:
    ch_versions = Channel.empty()
    ONCODRIVEFML(muts,    mutabs,    mutabind,    bedfile)

    emit:
    results = ONCODRIVEFML.out.tsv   // channel: [ val(meta), file(results) ]
    versions = ch_versions           // channel: [ versions.yml ]
}
