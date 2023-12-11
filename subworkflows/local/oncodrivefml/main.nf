include { ONCODRIVEFML } from '../../../modules/local/bbgtools/oncodrivefml/main'

workflow ONCODRIVEFML_ANALYSIS{

    take:
    muts
    mutabs
    mutabind
    bedfile

    main:

    ONCODRIVEFML(muts,    mutabs,    mutabind,    bedfile)

    emit:
    results = ONCODRIVEFML.out.tsv
    // plots = ONCODRIVEFML.out.plots

}
