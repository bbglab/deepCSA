process PLOT {

    tag "regressions"
    label 'process_single'

    container "docker.io/rblancomi/bbgregressions:dev"

    input:
    path config
    path data, stageAs: 'regressions/*'
    path metadata

    output:
    path ("plot/*")         , emit: plots
    path "versions.yml"     , topic: versions

    script:
    def predictors = params.bbgr_predictors.split(',').collect { it.trim() }
    """
    # update the YAML
    python3 <<EOF
    import yaml

    with open('${config}', 'r') as f:
        config_data = yaml.safe_load(f)

    config_data['general']['predictors_file'] = '${metadata.name}'
    config_data['general']['sample_column'] = '${params.bbgr_metadata_sampleIDcol}'
    config_data['general']['predictors'] = [${predictors.collect { "'$it'" }.join(', ')}]

    # Write modified config
    with open('updated_config.yaml', 'w') as f:
        yaml.dump(config_data, f, default_flow_style=False)
    EOF

    bbgregressions plot coefplot -config updated_config.yaml;

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
