include { SIGPROFILERASSIGNMENT } from '../../../modules/local/signatures/sigprofiler/assignment/main'


workflow SIGNATURES {
    take:
    // inputs
    matrix
    reference_signatures

    main:
    // actual code
    ch_versions = Channel.empty()
    SIGPROFILERASSIGNMENT(matrix, reference_signatures)

    // PROCESSDEPTHSTABLE()

    // PLOTDEPTHS()

    emit:
    depths   = SIGPROFILERASSIGNMENT.out.plots  // channel: [ val(meta), file(depths) ]
    versions = ch_versions                      // channel: [ versions.yml ]
}
