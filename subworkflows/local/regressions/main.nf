include { CREATE_INPUT          }    from '../../../modules/local/bbgtools/bbgregressions/create_input/main'
include { REGRESSIONS as MODELS }    from '../../../modules/local/bbgtools/bbgregressions/regressions/main'
include { PLOT                  }    from '../../../modules/local/bbgtools/bbgregressions/plot/main'


workflow REGRESSIONS{

    take:
    config
    data

    main:

    CREATE_INPUT(config, data)

    MODELS(config)

    PLOT(config)


    emit:
    inputs = CREATE_INPUT.out.inputs
    models = MODELS.out.models
    plots = PLOT.out.plots

}
