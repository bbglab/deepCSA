#!/usr/bin/env python3

"""
Run SigProfilerExtractor with the specified parameters.
"""

# Remember to use conda environment sigproext which has SigProfilerExtractor module installed

import click
from SigProfilerExtractor import sigpro as sig

@click.command()
@click.argument('input_type', default='matrix')  # When argument parameter is not needed to be specified in the command
@click.argument('output_dir')                     
@click.argument('input')                           
@click.option('--ref_genome', default='GRCh38', help='Reference genome to use') # When option parameter has to be specified in the command
@click.option('--max_sig', default=10, help='Maximum number of signatures')
@click.option('--nmf_replicates', default=100, help='Number of NMF replicates')
@click.option('--cpu', default=-1, help='Number of processors to be used to extract the signatures. Default value will use all available processors, which may cause a memory error.')

# Run SigProfilerExtractor
def main(input_type, output_dir, input, ref_genome, max_sig, nmf_replicates, cpu):
    
    sig.sigProfilerExtractor(
        input_type=input_type,
        output=output_dir,
        input_data=input, 
        reference_genome=ref_genome,
        minimum_signatures=1,
        maximum_signatures=max_sig,
        nmf_replicates=nmf_replicates,
	    cpu=cpu
    )

if __name__ == "__main__":
    main()

### EXAMPLE OF COMMAND: python3 sigprofiler_extractor2.py matrix ./output_dir count_matrix_wgs_20240821.txt --ref_genome GRCh38 --max_sig 10 --nmf_replicates 100 --cpu 10