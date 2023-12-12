include { ONCODRIVE3D } from '../../../modules/local/bbgtools/oncodrive3d/main'

workflow ONCODRIVE3D_ANALYSIS{

    take:
    muts
    mutabs
    mutabind

    main:

    ONCODRIVE3D(muts,    mutabs,    mutabind)

    emit:
    results     = ONCODRIVE3D.out.csv_genes
    results_pos = ONCODRIVE3D.out.csv_pos
    // plots = ONCODRIVE3D.out.plots

}
