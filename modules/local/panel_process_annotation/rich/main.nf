process POSTPROCESS_VEP_ANNOTATION {

    tag "${meta.id}"

    label 'cpu_low'
    label 'time_low'
    label 'process_high_memory'

    container "docker.io/ferriolcalvet/bgreference"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //         'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
    //         'biocontainers/pandas:1.5.2' }"


    input:
    tuple val(meta) , path(vep_annotated_file)

    output:
    tuple val(meta), path("*.compact.tsv") , emit: compact_panel_annotation
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assembly = task.ext.assembly ?: "hg38"
    def canonical_only = task.ext.canonical_only ? "True" : "False"
    // TODO
    // change panel postprocessing annotation into the same post processing annotation as before
    // keep it as the one for omega that is the one minimizing the computational processing
    // ensure that the columns are preserved whenever changing the version of VEP and
    // see if we want to separate information of the CANONICAL

    // FIXME
    // this selection of the columns from here looks very sensitive to changes in Ensembl, I think we should change it
    """
    zegrep -v '^##' ${vep_annotated_file} | cut -f 1,5,7,10,11,16,18,21-22 | awk '\$7!="-"' | uniq | \\
            tail -n +2 | \\
            awk -F'\\t' 'BEGIN {OFS = "\\t"} {split(\$1, a, "[_/]"); print a[1], a[2], a[3], a[4], \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9}' | \\
            gzip > ${prefix}.tmp.gz

    panel_postprocessing_annotation_rich.py \\
                    ${prefix}.tmp.gz \\
                    ${assembly} \\
                    ${vep_annotated_file.getBaseName()}.compact.tsv \\
                    ${canonical_only} ;
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${vep_annotated_file.getBaseName()}.compact.tsv;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
