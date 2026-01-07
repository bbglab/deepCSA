# Mutation catalog simulator

The main script

```
bash run_syntetic.sh
```

produces a collection of synthetic mutation catalogs (MAF and VCF format).

Each instance is specified by the following config parameters:

- Sample
- Set of genes
- dN/dS value
- Depth correction

Several random replicates are computed upon each config.

This information is provided in the config file (JSON format). For example, the following config file

```
{
    "samples": ["all_samples"],
    "omegas": [1.0, 2.0],
    "genes": ["ARID1A", "EGFR", "TP53"],
    "depth_corrections": [0.5, 1, 2],
    "n_replicates": 10
}
```

specifies a grid whereby 60 synthetic samples will be generated across 2 dN/dS and 3 depth correction values, with 10 random replicates in each case.

```depth_corrections``` allows the user to consider scenarios where the neutral mutability per site is rescaled, allowing simulation of different sparsity levels determined by lower/higher sequencing depths while keeping the same ground truth conditions. 

```samples``` represent reference samples in the deepCSA output from which the inferred neutral mutability per site is retrieved.