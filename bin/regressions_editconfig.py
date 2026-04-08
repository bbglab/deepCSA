#!/usr/bin/env python

import yaml
import json
import pandas as pd
import click

metric2filename = {"mutdensity": "all_mutdensities.tsv",
                    "omega": "all_omegas.tsv",
                    "omegagloballoc": "all_omegas_global_loc.tsv"}

@click.command()
@click.option('--config_file', type=click.Path(exists=True), required=True,
            help='Path to the config file to be edited')
@click.option('--mode', type=click.Choice(['custom', 'default']), required=True,
            help='deepCSA mode to run the regressions')
@click.option('--metric', type=click.Choice(['omega', 'omegagloballoc', "mutdensity"]), required=True,
            help='When mode=custom, specific metric to be included in this config file')
@click.option('--omega_res_file', type=click.Path(exists=True), required=True,
            help='When mode=default, omega results to define the genes for which to run regressions')
@click.option('--groups_file', type=click.Path(exists=True), required=True,
            help='When mode=default, groups information to only include single samples in the regressions')


def main(config_file: str,
        mode: str,
        metric: str,
        omega_res_file: str,
        groups_file: str):
    """
    Takes a bbgregressions formated config file and,
    depending on mode, updates it:
    *custom mode: retrieves the specific metric part of the file
    and saves again
    *default mode: updates elements and samples based on omega results
    and excluding group categories, respectively
    """

    # load config
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    # custom mode: config edited to contain metric info
    if mode == "custom":

        config_upd = {}
        config_upd["general"] = config["general"]
        config_upd["general"]["output_dir"] = "./"
        config_upd["plot"] = config["plot"]
        config_upd["metrics"] = {}
        print(metric)
        for config_metric in config["metrics"]:
            print(config_metric)
            print(config["metrics"][config_metric]["metric_name"])
            if config["metrics"][config_metric]["metric_name"] in metric: # this way to account for omega/omegaglolloc

                if metric == "mutdensity":
                    print("went through mutdensity")
                    config_upd["metrics"][config_metric] = config["metrics"][config_metric]
                    # file must correspond the name in the pipeline
                    config_upd["metrics"][config_metric]["file"] = metric2filename[metric]
                    with open(f'config_{mode}_{metric}.yaml', 'w') as f:
                        yaml.dump(config_upd, f, default_flow_style = False)
                    break

                elif metric == "omega":
                    print("went through omega")
                    if not config["metrics"][config_metric]["global_loc"]:
                        config_upd["metrics"][config_metric] = config["metrics"][config_metric]
                        # file must correspond the name in the pipeline
                        config_upd["metrics"][config_metric]["file"] = metric2filename[metric]
                        with open(f'config_{mode}_{metric}.yaml', 'w') as f:
                            yaml.dump(config_upd, f, default_flow_style = False)
                        break

                elif metric == "omegagloballoc":
                    print("went through omegaglobal")
                    print(config["metrics"][config_metric]["global_loc"])
                    if config["metrics"][config_metric]["global_loc"]:
                        print("when through global_loc = yes")
                        config_upd["metrics"][config_metric] = config["metrics"][config_metric]
                        # file must correspond the name in the pipeline
                        config_upd["metrics"][config_metric]["file"] = metric2filename[metric]
                        with open(f'config_{mode}_{metric}.yaml', 'w') as f:
                            yaml.dump(config_upd, f, default_flow_style = False)
                        break


    # default mode: elements are updated based on omega results; samples are updated to exclude groups
    elif mode == "default":

        # load omega res and select cohort-level significant unflagged genes (10 max)
        omega_res = pd.read_csv(omega_res_file, sep = "\t")
        omega_res = omega_res.loc[(omega_res["flagged"] == False)
                                & (omega_res["pvalue"] < 0.05)
                                & (omega_res["sample"] == "all_samples")
                                & (omega_res["impact"].isin(["truncating", "missense"]))
                                ]
        # selecting the
        omega_res_genes = omega_res.sort_values(by = "dnds", ascending = False).drop_duplicates(subset = "gene")["gene"].head(10)
        omega_res = omega_res[omega_res["gene"].isin(omega_res_genes)].reset_index(drop = True)
        if metric == "mutdensity":
            genes = list(omega_res["gene"].unique()) + ['ALL_GENES']
        else:
            omega_res["gene_impact"] = omega_res.apply(lambda row: f"{row['gene']}_{row['impact']}",
                                                        axis = 1)
            genes = list(omega_res["gene_impact"].unique())

        config["general"]["elements"] = genes

        # load 'all_samples' group and identify included samples
        with open(groups_file, 'r') as f:
            groups = json.load(f)
        samples = groups["all_samples"]
        config["general"]["samples"] = samples

        with open(f'config_{mode}_{metric}.yaml', 'w') as f:
            yaml.dump(config, f, default_flow_style = False)

if __name__ == '__main__':
    main()


