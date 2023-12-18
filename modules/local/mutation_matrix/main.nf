process COMPUTE_MATRIX {

    tag "$meta.id"
    label 'process_low'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container 'docker.io/ferriolcalvet/bgreference'

    input:
    tuple val(meta), path(mut_files)

    output:
    tuple val(meta), path("*.matrix.tsv")  , emit: matrix
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def filters = task.ext.filters ?: ""

    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    cat > mutations_subset.conf << EOF
    {
        ${filters}
    }
    EOF

    mut_profile.py matrix \\
                    --mut_file ${mut_files} \\
                    --out_matrix ${prefix}.matrix.tsv \\
                    --json_filters mutations_subset.conf \\
                    ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.matrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

// mutations_subset.conf \\

//     maf_df = maf_df.loc[
//     (maf_df["TYPE"] == "SNV")
//     & (maf_df["VAF"] <= 0.35)
//     & (maf_df["DEPTH"] >= 1000)
//     & (~maf_df["FILTER"].str.contains("not_in_panel"))            # only extended regions
//     & (~maf_df["FILTER"].str.contains("no_pileup_support"))
//     & (~maf_df["FILTER"].str.contains("low_complex_repetitive"))  # discard also these regions
//     & (~maf_df["FILTER"].str.contains("low_mappability"))         # discard also these regions
//     & (~maf_df["FILTER"].str.contains("n_rich"))
//     & (~maf_df["FILTER"].str.contains("repetitive_variant"))
//     & (~maf_df["FILTER"].str.contains("other_sample_SNP"))
// ]
