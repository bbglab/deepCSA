#dry run
#snakemake  --use-conda --snakefile create_signatures_extraction_input.snk.py -np
#run
#snakemake  --use-conda --snakefile create_signatures_extraction_input.snk.py

#deepCSA_run = "/data/bbg/nobackup/prominent/kidney/analysis/2024-11-18_LCMs"
#OUTPUT_FOLDER = f"/data/bbg/nobackup/prominent/kidney/analysis/2024-11-18_LCMs/input_for_signature_extraction"
#deepCSA_run = "/data/bbg/nobackup/prominent/kidney/analysis/2024-09-30_deepCSA_LCMs"
#OUTPUT_FOLDER = f"/data/bbg/nobackup/prominent/kidney/analysis/2024-09-30_deepCSA_LCMs/input_for_signature_extraction"
deepCSA_run = "/data/bbg/nobackup/lung_duplex/analysis/fullcohortnormal/2024-12-30_run_all"
OUTPUT_FOLDER = f"{deepCSA_run}/input_for_signature_extraction"


SUB_DIR_1 = f"{OUTPUT_FOLDER}/results"
SUB_DIR_2 = f"{OUTPUT_FOLDER}/logs"
COMMON_CONDA_ENV="hdp_env"

SCRIPTS_PATH = "./HDP_input/bin"
import subprocess, sys, os, glob

# Input parameters  ---------------------------------

# Rules             ---------------------------------

# Ensure the output folder exists

rule all:
    input:
        f"{SUB_DIR_1}/Mutation_stats_per_genome_per_sample.med.all.csv",
        f"{SUB_DIR_1}/Table_for_SigProfilerExtractor.med.all.csv",
        f"{SUB_DIR_1}/count_matrix_hdp.med.all.csv",
        f"{SUB_DIR_1}/count_matrix_hdp.med.all.rds",
        f"{SUB_DIR_1}/treeLayer_hdp.med.all.csv",
        f"{SUB_DIR_1}/treeLayer_hdp.med.all.rds"

rule COUNT_MUTATIONS:
    input:
        input_profile=f"{deepCSA_run}/computeprofile/",
        input_mutations=f"{deepCSA_run}/somaticmutations/",
        input_samples_json=f"{deepCSA_run}/table2group/samples.json"
    output:
        f"{SUB_DIR_1}/Mutation_stats_per_genome_per_sample.med.all.csv"
    params:
        create_dir=OUTPUT_FOLDER
    conda:
        COMMON_CONDA_ENV  # Specify the correct path to your environment
    log:
        f"{SUB_DIR_2}/count_mutations.log"  # Specify a log file
    shell:
        """
        Rscript --vanilla {SCRIPTS_PATH}/count_mutations.R {input.input_profile} {input.input_mutations} {output} {input.input_samples_json} > {log} 2>&1
        """

rule SIGPROFILER_TABLE:
    input:
        input_profile=f"{deepCSA_run}/computeprofile/",
        input_samples_json=f"{deepCSA_run}/table2group/samples.json"
    output:
        f"{SUB_DIR_1}/Table_for_SigProfilerExtractor.med.all.csv"
    params:
        create_dir=OUTPUT_FOLDER
    conda:
        COMMON_CONDA_ENV  # Specify the correct path to your environment YAML
    log:
        f"{SUB_DIR_2}/create_sigprofiler_table.log"  # Specify a log file
    shell:
        """
        Rscript --vanilla {SCRIPTS_PATH}/create_table_for_SigProfilerExtractor.R {input.input_profile} {output} {input.input_samples_json} > {log} 2>&1
        """


rule HDP_INPUT:
    input:
        mutatation_counts=f"{SUB_DIR_1}/Table_for_SigProfilerExtractor.med.all.csv"
    output:
        out_csv=f"{SUB_DIR_1}/count_matrix_hdp.med.all.csv",
        out_rds=f"{SUB_DIR_1}/count_matrix_hdp.med.all.rds"
    params:
        create_dir=OUTPUT_FOLDER
    conda:
        COMMON_CONDA_ENV  # Specify the correct path to your environment YAML
    log:
        f"{SUB_DIR_2}/create_hdp_input.log"  # Specify a log file
    shell:
        """
        Rscript --vanilla {SCRIPTS_PATH}/create_input_for_hdp.R {input.mutatation_counts} {output.out_csv} {output.out_rds} > {log} 2>&1
        """

rule HDP_TREE:
    input:
        f"{SUB_DIR_1}/Mutation_stats_per_genome_per_sample.med.all.csv"
    output:
        out_csv=f"{SUB_DIR_1}/treeLayer_hdp.med.all.csv",
        out_rds=f"{SUB_DIR_1}/treeLayer_hdp.med.all.rds"
    params:
        create_dir=OUTPUT_FOLDER
    conda:
        COMMON_CONDA_ENV  # Specify the correct path to your environment YAML
    log:
        f"{SUB_DIR_2}/create_treeLayer_for_hdp.log"  # Specify a log file
    shell:
        """
        Rscript --vanilla {SCRIPTS_PATH}/create_treeLayer_for_hdp.R {input} {output.out_csv} {output.out_rds} > {log} 2>&1
        """

rule create_output_directory:
    output:
        directory(OUTPUT_FOLDER)  # Declare OUTPUT_FOLDER as a directory output
    shell:
        """
        mkdir -p {output} {SUB_DIR_1} {SUB_DIR_2}
        """


