
include { PLOT_VAF_MUTDENSITY_DEPTH    as PLOTVAFMUTDENSITY     } from '../../../modules/local/plot/vaf_mutdensity_depth/main'
include { PLOT_SELECTION_DEPTH          as PLOTSELECTIONDEPTH    } from '../../../modules/local/plot/selection_depth/main'


workflow PLOT_DEPTH_RELATIONSHIPS {

    take:
    somatic_mutations       // MAF files with mutations
    all_mutdensities        // Mutation density files
    depth_per_gene          // Depth per gene files  
    omega_results           // Omega results (optional)
    oncodrivefml_results    // OncodriveFML results (optional)


    main:

    // Prepare channels for VAF and mutation density plotting
    // Combine MAF, mutation density, and depth files
    somatic_mutations
    .map { meta, maf -> tuple(meta.id, meta, maf) }
    .set { maf_ch }
    
    all_mutdensities
    .map { meta, mutdens -> tuple(meta.id, mutdens) }
    .set { mutdens_ch }
    
    depth_per_gene
    .map { meta, depth -> tuple(meta.id, depth) }
    .set { depth_ch }
    
    // Join all inputs for VAF/mutation density plotting
    maf_ch
    .join( mutdens_ch, remainder: true )
    .join( depth_ch, remainder: true )
    .map { id, meta, maf, mutdens, depth ->
        // Handle missing files by providing NO_FILE placeholder
        def maf_file = maf ?: file('NO_FILE')
        def mutdens_file = mutdens ?: file('NO_FILE')
        def depth_file = depth ?: file('NO_FILE')
        tuple(meta, maf_file, mutdens_file, depth_file)
    }
    .set { vaf_mutdens_input }
    
    PLOTVAFMUTDENSITY(vaf_mutdens_input)
    
    // Prepare channels for selection metric plotting if available
    if ( omega_results && oncodrivefml_results ) {
        omega_results
        .map { meta, omega -> tuple(meta.id, omega) }
        .set { omega_ch }
        
        oncodrivefml_results
        .map { meta, ofml -> tuple(meta.id, ofml) }
        .set { ofml_ch }
        
        // Join omega, oncodrivefml, and depth
        omega_ch
        .join( ofml_ch, remainder: true )
        .join( depth_ch, remainder: true )
        .map { id, omega, ofml, depth ->
            // Get meta from one of the channels
            def meta = [id: id]
            def omega_file = omega ?: file('NO_FILE')
            def ofml_file = ofml ?: file('NO_FILE')
            def depth_file = depth ?: file('NO_FILE')
            tuple(meta, omega_file, ofml_file, depth_file)
        }
        .filter { meta, omega, ofml, depth -> 
            !depth.name.equals('NO_FILE') && (!omega.name.equals('NO_FILE') || !ofml.name.equals('NO_FILE'))
        }
        .set { selection_input }
        
        PLOTSELECTIONDEPTH(selection_input)
    }

    emit:
    vaf_mutdensity_plots    = PLOTVAFMUTDENSITY.out.plots
    selection_plots         = omega_results && oncodrivefml_results ? PLOTSELECTIONDEPTH.out.plots : Channel.empty()

}
