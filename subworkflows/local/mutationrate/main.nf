include { MUTRATE as MUTRATE_ALL} from '../../../modules/local/computemutrate/main'
include { MUTRATE as MUTRATE_PROTAFFECT} from '../../../modules/local/computemutrate/main'
include { MUTRATE as MUTRATE_NONPROTAFFECT} from '../../../modules/local/computemutrate/main'


workflow MUTATION_RATE{
    take:
    maf
    depths
    consensus_annot_panel_all
    consensus_annot_panel_protaffect
    consensus_annot_panel_nonprotaffect

    main:
    ch_versions = Channel.empty()

    MUTRATE_ALL(maf, depths, consensus_annot_panel_all)
    MUTRATE_PROTAFFECT(maf, depths, consensus_annot_panel_protaffect)
    MUTRATE_NONPROTAFFECT(maf, depths, consensus_annot_panel_nonprotaffect)

    // PLOTMUTRATE()

    emit:
    mutrates_all   = MUTRATE_ALL.out.mutrates
    mutrates_protaffect   = MUTRATE_PROTAFFECT.out.mutrates
    mutrates_nonprotaffect   = MUTRATE_NONPROTAFFECT.out.mutrates

    versions = ch_versions                // channel: [ versions.yml ]
}
