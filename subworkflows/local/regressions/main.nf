include { RUNREGRESSIONS    as RUNREGRESSIONS             } from '../../../modules/local/runregressions/main'

workflow REGRESSIONS{
    take:
    features_table

    main:
    ch_versions = Channel.empty()

    RUNREGRESSIONS(features_table)
    ch_versions = ch_versions.mix(RUNREGRESSIONS.out.versions)

    emit:
    regressions = RUNREGRESSIONS.out.regressions
    versions = ch_versions                // channel: [ versions.yml ]
}
