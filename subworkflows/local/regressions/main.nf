include { EDITCONFIG            }    from '../../../modules/local/bbgtools/bbgregressions/editconfig/main'
include { CREATE_INPUT          }    from '../../../modules/local/bbgtools/bbgregressions/create_input/main'
include { REGRESSIONS as MODELS }    from '../../../modules/local/bbgtools/bbgregressions/regressions/main'
include { PLOT                  }    from '../../../modules/local/bbgtools/bbgregressions/plot/main'

workflow REGRESSIONS{

    take:
    config
    data
    metadata
    mode
    metric
    omega_res
    groups

    main:

    EDITCONFIG(config, mode, metric, omega_res, groups)

    CREATE_INPUT(EDITCONFIG.out.config, data)

    MODELS(EDITCONFIG.out.config, CREATE_INPUT.out.inputs, metadata)

    PLOT(EDITCONFIG.out.config, MODELS.out.models, metadata)


    emit:
    inputs = CREATE_INPUT.out.inputs
    models = MODELS.out.models
    plots = PLOT.out.plots

}
