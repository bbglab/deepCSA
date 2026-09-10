process COMPUTE_PROFILE {

    tag "$meta.id"
    label 'cpu_low'
    label 'mem_low'


    label 'deepcsa_core'

    input:
    tuple val(meta) , path(matrix), path(trinucleotide)
    path( wgs_trinucleotides )
    tuple val(meta2), path( cohort_profile , stageAs: 'global_mutprofile.tsv')

    output:
    tuple val(meta), path("*.profile.tsv")                                  , emit: profile
    tuple val(meta), path("*.proportion_mutations.tsv")     , optional:true , emit: panel_proportions
    tuple val(meta), path("*.proportion_mutations.WGS.tsv") , optional:true , emit: wgs_proportions
    tuple val(meta), path("*.matrix.WGS.tsv")               , optional:true , emit: wgs
    tuple val(meta), path("*.matrix.WGS.sigprofiler.tsv")   , optional:true , emit: wgs_sigprofiler
    tuple val(meta), path("*.profile_stability.tsv")                        , emit: profile_stability

    tuple val(meta), path("*.pdf")                          , optional:true , emit: plots

    path "versions.yml"                                                     , topic: versions


    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def wgs_trinuc = wgs_trinucleotides ? "--wgs --wgs_trinucleotide_counts ${wgs_trinucleotides}" : ""
    def smoothing_args = task.ext.smoothing ? "--smoothed --prior_profile ${cohort_profile}" : ""
    """
    mut_profile.py profile \\
                    --sample_name ${prefix} \\
                    --mutation_matrix ${matrix} \\
                    --trinucleotide_counts ${trinucleotide} \\
                    ${wgs_trinuc} \\
                    ${smoothing_args} \\
                    ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """


    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.profile.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
