process REGRESSIONS {

    tag "regressions"
    label 'process_single'

    container "docker.io/bbglab/bbgregressions:0.1.0"

    input:
    path config

    output:
    path ("regressions/*")       , emit: models
    path "versions.yml"          , topic: versions

    script:
    def predictors = params.bbgr_predictors.split(',').collect { it.trim() }
    """
    # update the YAML
    python3 <<EOF
    import yaml

    with open('${config}', 'r') as f:
        config_data = yaml.safe_load(f)

    config_data['predictors_file'] = '${params.bbgr_metadata}'
    config_data['sample_column'] = '${params.bbgr_metadata_sampleIDcol}'
    config_data['predictors'] = [${predictors.collect { "'$it'" }.join(', ')}]

    # Write modified config
    with open('${config}', 'w') as f:
        yaml.dump(config_data, f, default_flow_style=False)
    EOF

    bbgregressions regressions -config modified_config.yml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}


