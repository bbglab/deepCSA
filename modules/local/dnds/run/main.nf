process RUN_DNDS {
    tag "$meta.id"


    container 'docker.io/ferriolcalvet/dnds:latest'

    input:
    tuple val(meta) , path(mutations_table), path(depths)
    tuple val(meta2), path(ref_cds)
    path (covariates)

    output:
    tuple val(meta), path("*.cv.tsv*")          , emit: results_cv
    tuple val(meta), path("*.loc.tsv*")         , emit: results_local
    tuple val(meta), path("*.globaldnds.tsv*")  , emit: results_global    
    path "versions.yml"                         , topic: versions


    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    dNdS_run.R --inputfile ${mutations_table} \\
                --outputprefix ${prefix} \\
                --samplename ${prefix} \\
                --covariates ${covariates} \\
                --referencetranscripts ${ref_cds} \\
                --genedepth ${depths} \\
                ${args}
    # --cores ${task.cpus}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dNdScv: 0.1.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.cv.tsv
    touch ${prefix}.loc.tsv
    touch ${prefix}.globaldnds.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dNdScv: 0.1.0
    END_VERSIONS
    """
}


// --snvsonly"), type="logical",
//             default=FALSE,
//             help="Only use SNVs for the analysis [default= %default]", metavar="logical")
