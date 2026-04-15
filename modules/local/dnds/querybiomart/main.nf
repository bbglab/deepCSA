process QUERY_BIOMART {

    tag "$meta.id"
    label 'cpu_single_fixed'
    label 'time_low'
    label 'process_high_memory'

    container 'docker.io/ferriolcalvet/querybiomart:0.2.0'

    input:
    tuple val(meta) , path(panel)

    output:
    tuple val(meta), path("cds_biomart.txt"), emit: complete_biomart
    path "versions.yml"                                 , topic: versions

    script:
    """
    cut -f 6 ${panel} | tail -n +2 | sort -u | awk '\$1!="-"' | tr -s '\\n' ',' | sed 's/,\$//' > genes_list.txt

    biomart_query.R --genelist genes_list.txt \\
                            --outputfile cds_biomart.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Ensembl BioMart: v111
    END_VERSIONS
    """


    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch cds_biomart.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Ensembl BioMart: v111
    END_VERSIONS
    """
}

// FIXME: right now the ensembl archive version of biomart is hardcoded
// we should define some map and make sure that it is updated accordingly