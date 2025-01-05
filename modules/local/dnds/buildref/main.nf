process BUILD_REFCDS {

    tag "$meta.id"
    label 'cpu_single_fixed'
    label 'time_low'
    label 'process_high_memory'


    container 'docker.io/ferriolcalvet/dnds:latest'

    input:
    tuple val(meta) , path(biomart_cds)
    tuple val(meta2), path(reference_genome)
    

    output:
    tuple val(meta), path("RefCDS.rda") , emit: ref_cds
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Rscript "library(dndscv); buildref(\"${biomart_cds}\", \"${reference_genome}\", outfile = \"RefCDS.rda\")"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R     : 
        dNdScv: 0.1.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch RefCDS.rda

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R     : 
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

// process runRScript {
//     input:
//     path 'input.txt'

//     output:
//     path 'output.txt'

//     script:
//     """
//     #!/usr/bin/env Rscript

//     # Your R code here
//     data <- read.table("input.txt")
//     # Process data
//     write.table(processed_data, "output.txt")
//     """
// }