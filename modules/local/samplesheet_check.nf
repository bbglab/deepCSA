process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in bbglab/deepCSA/bin/
    def s3_params = params.mountS3 ? "--s3startingPoint ${params.s3startingPoint} --s3bucketName ${params.s3bucketName}" : ""


    """
    check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv \\
        $s3_param

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
