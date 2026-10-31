#!/usr/bin/env python

import click
import pandas as pd

from read_utils import custom_na_values


@click.command()
@click.option("--all-samples-file", required=True, type=click.Path(exists=True))
@click.option("--mutations-file", required=True, type=click.Path(exists=True))
@click.option("--output-file", required=True, type=click.Path())
@click.option("--percentile", default=0.995, type=click.FloatRange(0.0, 1.0))
def main(all_samples_file, mutations_file, output_file, percentile):
    """Add cohort-capped ALT_DEPTH as EFFECTIVE_DEPTH to a mutations TSV."""
    all_samples = pd.read_csv(
        all_samples_file, sep="\t", header=0, na_values=custom_na_values
    )
    mutations = pd.read_csv(
        mutations_file, sep="\t", header=0, na_values=custom_na_values
    )

    alt_depth = pd.to_numeric(all_samples["ALT_DEPTH"], errors="raise")
    cap = alt_depth.quantile(percentile)
    if pd.isna(cap):
        raise click.ClickException(
            f"{all_samples_file} has no numeric ALT_DEPTH values to calculate a cap"
        )

    mutations["EFFECTIVE_ALT_DEPTH"] = pd.to_numeric(
        mutations["ALT_DEPTH"], errors="raise"
    ).clip(upper=cap)
    mutations.to_csv(output_file, sep="\t", header=True, index=False)
    click.echo(
        f"Added EFFECTIVE_ALT_DEPTH using the {percentile:.3%} ALT_DEPTH percentile "
        f"({cap:g}) to {output_file}"
    )


if __name__ == "__main__":
    main()
