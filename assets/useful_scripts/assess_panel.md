# Experimental design principles

## Definition of relationship between desired mutations and required depth

Given a fixed panel size based on which are the regions of interest, a desired number of mutations (enough to get good selection/mutagenesis metrics) and a known estimation of the expected mutation density of the sample, you can estimate which is the sequencing depth required using this formula:

$$
\text{Required depth} = \frac{\text{minimum desired mutations}}{(\text{region of interest (bp)}) \times (\text{expected mutation density (muts/bp)})}
$$

To get an estimation of the expected mutation density you can use this formula even if it is with some approximate values:

$$
\text{expected mutation density (muts/bp)} = \text{mutation rate (mutations} \times \text{genome} \times \text{year)} \times \text{time (year)}
$$

For the mutation rate this could be a possible formula; but you may also rely on previously published data.

$$
\text{mutation rate (} \frac{\text{mutations}}{\text{bp * year}} \text{)}  = \frac{\text{observed mutations}}{ \text{genome} \times \text{age (years)} }
$$

In this example, we are assuming that the units of the mutation rate are mutations x genome x year, but this could be in any other units
Particularly the time units could be different from year.

In case that the sequencing depth is restricted by the availability of DNA or any other reason, you can decide to adjust the panel so that by increasing the depth you can still reach a stable measurement of the mutation density and other variables.

An additional scenario to consider is the trinucleotide mutation probabilities. In a given sample, the different sites may have a very different mutation probablility depending on which are the active mutational processes. If these are known, you can factor in this information to adjust the requirements of the experimental design.
