process FILTERBED {

    tag "$meta.id"
    label 'process_high'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta)     , path(maf)
    tuple val(meta2)    , path(bedfile)

    output:
    tuple val(meta), path("*.tsv.gz")  , emit: maf
    path "versions.yml"                , topic: versions


    script:
    def filtername = task.ext.filtername ?: "covered"
    def positive = task.ext.positive ?: false
    def positive_flag = positive ? "--positive" : ""

    """
    filterbed.py \\
        --sample-maf-file ${maf} \\
        --bedfile ${bedfile} \\
        --filtername ${filtername} \\
        ${positive_flag};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.vep.summary.tab.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}