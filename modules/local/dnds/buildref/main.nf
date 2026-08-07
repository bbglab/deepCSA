process BUILD_REFCDS {

    tag "$meta.id"


    container 'docker.io/ferriolcalvet/dnds:latest'

    input:
    tuple val(meta) , path(biomart_cds)
    path(reference_genome)


    output:
    tuple val(meta), path("RefCDS_custom.rda") , emit: ref_cds
    path "versions.yml"                        , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Rscript -e "library(dndscv); buildref('${biomart_cds}', '${reference_genome}', outfile = 'RefCDS_custom.rda')"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R     : 4.4.2
        dNdScv: 0.1.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch RefCDS_custom.rda

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R     : 4.4.2
        dNdScv: 0.1.0
    END_VERSIONS
    """
}