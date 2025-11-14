# bbglab/deepCSA: Tools explanation

## Table of contents

- [Introduction](#introduction)

## Introduction

Here, you can find an explanation of the different computations, tools or metrics implemented in deepCSA.

## Publications with detailed explanation

We are in the process of completing the documentation, but in the meantime you can check the recently published [paper and its supplementary material for more details](https://www.nature.com/articles/s41586-025-09521-x).

## Adjusted mutation density

This can be found in the `mutdensityadj` folder of the deepCSA output.

### Goal

The goal is to define a mean mutation density estimate for a specific subset of mutation sites. For example, we might be interested in measuring what is the mutation density at the set of missense mutation sites -- i.e., mutation sites in the CDS that induce missense mutations.

### Motivation

In the context of bulk-ultradeep sequencing, there are several factors that can influence the occurrence of mutations at a given site.

- trinucleotide context of the site

- neutral mutagenesis, defined as:

  - normalized profile: vector of relative mutabilities for each trinucleotide context

  - exposure associated with the profile, expressed as mutation burden

- sequencing depth of the site

- selection of the site

We would like to define a way to compute the mutation density that corrects for potential confounders like depth and triplet content, thus rendering the estimates more comparable. For example, the missense sites in gene A may imply that more missense mutations per site are observed than in gene B, even if both genes are subject to the same exact mutational process.

Given a collection $\mathcal{C}=\\{i\\}$ of mutation sites (e.g. missense mutations) and a number of mutations $m$ observed in this context, we would like to come up with an expression of the type $m/L$, where $L$ denotes the effective number of sites. 

Assuming even depth $D$ for simplicity, the current *flat method* would simply take :

```math
L=\frac{1}{3}\cdot|\mathcal{C}|\cdot D
```

### Assumptions

Throughout we make the following assumptions:

- the only way mutations come about is by neutral mutagenesis

- discrepancies between neutral mutagenesis and the observed mutations are due to selection

### Definitions/Notation

Mutation sites are defined as specific single base nucleotide changes. Each mutation site corresponds to a tri-nucleotide context of the form $F_5RF_3>A$ where $F_5, F_3$ are the 5' and 3' flanking nucleotides, respectively, R is the reference nucleotide and A is the alternate allele. Therefore, there are 3 mutation sites per genomic position.

Let $\\{P_c\\}$ be the relative mutabilities per trinucleotide context $c$: given a mutation site with tri-nucleotide context $c_1$, we expect to observe $P_{c_1}/P_{c_2}$ more mutations than in a mutation site with triplet context $c_2$. For a more univocal definition, we may want to impose the following condition:

```math
\sum_c P_c = 1
```

Note that in order to render absolute mutabilities, the $P_c$ should be multiplied by a suitable scaling factor $\alpha$.

Let $D_i$ be the sequencing depth at mutation site $i$.

Let $m$ denote the total number of mutations observed in the target class -- i.e., mapping to the prescribed collection of mutation sites where we are interested to compute the mutation density.

For each mutation site $i$, we will denote $c(i)$ its triplet context.

### Method description

Given a collection $\mathcal{C}=\\{i\\}$ of mutation sites (e.g. missense mutations) and a number of mutations $m$ observed in this context, we want to compute an effective length $L$ that represents the entire collection of sites

```math
L=\sum_{i\in\mathcal{C}} P_{c(i)}\cdot D_i
```

Note that this definition of length has an arbitrary interpretation, because of the multiple possible choices that we can make for the relative mutational profile. To render a universal length, we must set a scale to represent the mutational profile with a concrete meaning. What would be a good practical reference?

One possibility would be to choose a scaling factor $\alpha$ that represents an absolute mutability equivalent to an expected 1 mutation per megabase in one mappable copy of the genome.

If $n_c$ is the number of mutation sites with trinucleotide context $c$ comprised in the mappable genome and $G$ is the total length of the mappable genome (expressed in nucleotides), $\alpha$ must be such that:

```math
\alpha\sum_c n_cP_c = G / 10^6
```

Because $G=\frac{1}{3}\sum_c n_c$, we can compute $\alpha$:

```math
\hat\alpha = \frac{G}{10^6\cdot\sum_c n_cP_c} = \frac{1}{3\cdot 10^6} \frac{\sum_c n_c}{\sum_c n_cP_c}
```

Applying this scaling, the effective length would take the following form:

```math
L=\hat\alpha\sum_i P_{c(i)}\cdot D_i = \hat\alpha\sum_c \left( P_c\cdot \sum_{i\in c} D_i \right)
```

### Units

Defining the mutation density as $m/L$, according to the preceding explanation, the units in which the readout is expressed are the following: mutations / Mb sequenced.

## Omega

For more explanations on omega go to the [corresponding repo](https://github.com/bbglab/omega).

## Site comparison

The site comparison step takes advantage of the computation of mutabilities in [omega](https://github.com/bbglab/omega), and then compares these mutabilities either by residue, residue change or nucleotide change.

## Mutational signatures

We provide two different strategies for signature analysis.

- Using SigProfilerAssignment with a set of known SBS signatures

- Using a Hierarchical Dirichlet Process algorithm developed by Nicola Robets and compacted by the McGranahan lab into a wrapped version.

Additionally one could run SigProfilerExtractor on the data but this needs to be done externally.
