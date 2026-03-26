# Mutation catalog simulator

The main script

```
bash run_syntetic.sh
```

produces a collection of synthetic mutation catalogs (MAF and VCF format).

Each instance is specified by the following config parameters:

- Sample
- Gene
- dN/dS value
- Consequence type
- Depth correction

Several random replicates are computed upon each config.

This information is provided in the config file (JSON format). For example, the following config file

```
{
    "samples": ["all_samples"],
    "omegas": [1.0, 2.0],
    "csqn_types": ["missense", "nonsense", "nonsynonymous"]
    "genes": ["ARID1A", "EGFR", "TP53"],
    "depth_corrections": [0.5, 1, 2],
    "n_replicates": 10
}
```

specifies a grid whereby 120 synthetic samples will be generated across for two dN/dS values, three consequence types, and three depth correction values, with 10 random replicates in each case.

```depth_corrections``` is a scaling factor that allows the user to consider scenarios where the neutral mutability per site is rescaled, allowing simulation of different sparsity levels determined by lower/higher sequencing depths while keeping the same ground truth conditions.

The typical usecase is to simulate samples that are representative of the average sample in terms of depth. In that case,
the recommended depth correction is simply the inverse of the number of samples if the reference sample taken in "all_samples".

```samples``` represent reference samples in the deepCSA output from which the inferred neutral mutability per site is retrieved.