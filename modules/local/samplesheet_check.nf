process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    path samplesheet
    val require_bams

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", topic: versions


    script: // This script is bundled with the pipeline, in bbglab/deepCSA/bin/
    def validate_bams = require_bams ? '--bam-required' : ''
    """
    check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv \\
        ${validate_bams}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
