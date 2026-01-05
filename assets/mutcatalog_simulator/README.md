# Mutation catalog simulator

The main script

```
bash run_syntetic.sh
```

produces a collection of synthetic mutation catalogs (MAF and VCF format).

Each instance is specified by the following parameters:

- Sample
- Gene
- dN/dS
- Consequence type ("nonsense", "missense")
- Depth correction
- Random replicate

This information is provided in the config file (JSON format). For example, the following config file

```
{
    "samples": ["all_samples"],
    "omegas": [1.0, 2.0],
    "genes": ["ARID1A", "EGFR", "TP53"],
    "csqn_types": ["missense","nonsense"],
    "depth_corrections": [1],
    "n_replicates": 10
}
```

specifies a grid whereby 120 synthetic samples will be generated across two dN/dS, three genes and two consequence types, with 10 random replicates in each case.

```depth_corrections``` allows the user to consider scenarios where the average depth of the sequencing is lower $d<1$ or higher $d>1$ than the actual depth reported in ```samples```.