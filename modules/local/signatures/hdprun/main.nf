process RUN_HDP_WRAPPER {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/spellegrini87/hdp:1.0.0'

    input:
    tuple val(meta) , path(matrix)
    tuple val(meta2), path(treelayers)
    path(ref_signatures)

    output:
    tuple val(meta), path("**.pdf")     , emit: plots
    tuple val(meta), path("**.csv")     , emit: stats
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
     # First, create an R script
    cat <<EOF > process_data.R
    data = read.table("${matrix}", header = FALSE)
    rownames(data) <- data[,c(1)]
    colnames(data) <- data[c(1),]
    data <- data[-c(1),]
    data[,-c(1)] <- sapply(data[,-c(1)], as.numeric)
    data <- data[,-c(1)]
    saveRDS(data, file = "${prefix}.hdp.rds")
    EOF

    # Run the R script
    Rscript process_data.R


    cat <<-END_INPUT > input.txt
    ${prefix}.hdp.rds\t${treelayers}
    END_INPUT
    
    HDP_wrapper.sh \\
                input.txt \\
                ${ref_signatures} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R       : \$(Rscript --version | sed -e 's/.*version //g')
        Rscript : \$(Rscript --version | sed -e 's/.*version //g')
        mSigHdp : 2.1.2
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R       : \$(Rscript --version | sed -e 's/.*version //g')
        Rscript : \$(Rscript --version | sed -e 's/.*version //g')
        mSigHdp : 2.1.2
    END_VERSIONS
    """
}
