process MERGE_VEP {
    label 'process_medium'

    conda "conda-forge::coreutils=8.31"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'https://depot.galaxyproject.org/singularity/coreutils:8.31--h14c3975_0' : 'quay.io/biocontainers/coreutils:8.31--h14c3975_0' }"

    input:
    tuple val(meta), path(vep_files)

    output:
    tuple val(meta), path("*.merged.tab.gz"), emit: tab
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Output filename
    out_file="${prefix}.merged.tab"

    # Sort input files by chromosome (e.g. chr1, chr2, ..., chrX, chrY)
    sorted_files=\$(printf "%s\\n" ${vep_files} | sort -V)

    # Initialize header_written flag
    header_written=false

    # Merge files
    for f in \$sorted_files; do
        echo "Processing \$f..."
        if [ "\$header_written" = false ]; then
            zcat "\$f" | awk '{if(\$0 ~ /^##/) print; else if(\$0 ~ /^#/) {print; next} else print}' > "\$out_file"
            header_written=true
        else
            zcat "\$f" | awk 'NR>1 {print}' >> "\$out_file"
        fi
    done

    # Compress final merged output
    gzip "\$out_file"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(sort --version | head -n1 | sed 's/.*) //; s/ .*\$//')
        gzip: \$(gzip --version | head -n1 | sed 's/.*) //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.tab.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(sort --version | head -n1 | sed 's/.*) //; s/ .*\$//')
        gzip: \$(gzip --version | head -n1 | sed 's/.*) //; s/ .*\$//')
    END_VERSIONS
    """
}