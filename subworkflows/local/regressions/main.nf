include { RUNREGRESSIONS    as RUNREGRESSIONS             } from '../../../modules/local/runregressions/main'

workflow REGRESSIONS{
    take:
    regressions_config

    main:
    ch_versions = Channel.empty()

    RUNREGRESSIONS(regressions_config)
    ch_versions = ch_versions.mix(RUNREGRESSIONS.out.versions)

    emit:
    regressions = RUNREGRESSIONS.out.regressions
    versions = ch_versions                // channel: [ versions.yml ]
}
