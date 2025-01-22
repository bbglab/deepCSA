process RUN_DNDS {
    tag "$meta.id"
    label 'cpu_single_fixed'
    label 'time_low'
    label 'process_high_memory'


    container 'docker.io/ferriolcalvet/dnds:latest'

    input:
    tuple val(meta) , path(mutations_table), path(depths)
    tuple val(meta2), path(ref_cds)
    path (covariates)

    output:
    tuple val(meta), path("*.out.tsv*") , emit: results
    path "versions.yml"                 , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dNdS_run.R --inputfile ${mutations_table} \\
                --outputfile ${prefix}.out.tsv \\
                --samplename ${prefix} \\
                --covariates ${covariates} \\
                --referencetranscripts ${ref_cds} \\
                --genedepth ${depths}
    # --cores ${task.cpus}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dNdScv: 0.1.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.out.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dNdScv: 0.1.0
    END_VERSIONS
    """
}

// "--referencetranscripts"
// default="/workspace/projects/prominent/analysis/dNdScv/data/reference_files/RefCDS_human_latest_intogen.rda",
// --covariates
//     "/workspace/projects/prominent/analysis/dNdScv/data/reference_files/covariates_hg19_hg38_epigenome_pcawg.rda",
//             help="Human GRCh38 covariates file [default= %default]", metavar="character"),
// --genelist"), type="character",
//             default=NULL,
//             help="Gene list file [default= %default]", metavar="character"),
// --genedepth"), type="character",
//             default=NULL,
//             help="Gene depth file (2 columns: GENE\tAVG_DEPTH) [default= %default]", metavar="character"),
// --snvsonly"), type="logical",
//             default=FALSE,
//             help="Only use SNVs for the analysis [default= %default]", metavar="logical")
