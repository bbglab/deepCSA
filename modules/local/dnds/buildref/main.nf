process BUILD_REFCDS {

    tag "$meta.id"
    label 'cpu_single_fixed'
    label 'time_low'
    label 'process_high_memory'


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
    """
    Rscript -e "library(dndscv); buildref('${biomart_cds}', '${reference_genome}', outfile = 'RefCDS_custom.rda')"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R     : 4.4.2
        dNdScv: 0.1.0
    END_VERSIONS
    """

    stub:
    """
    touch RefCDS_custom.rda

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R     : 4.4.2
        dNdScv: 0.1.0
    END_VERSIONS
    """
}