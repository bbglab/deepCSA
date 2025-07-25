#!/opt/conda/bin/Rscript --vanilla

# standalone version of deepCSA depth script, to be used in the deepCSA pipeline
# dndscv-functions dndsloc and fit_substmodel taken from dNdScv package: https://github.com/im3sanger/dndscv
library(data.table)
library(tidyverse)
library(abind)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(IRanges)
options(dplyr.summarise.inform = FALSE)

# load dndscv substmodel
substmodel = read.delim("/dndscv_table/submod_192r_3w.tsv") |> as.matrix()

args = commandArgs(trailingOnly = TRUE)
# define output list elements
output_list = plot_list = list()

# Uncomment following lines for testing/checking
# data("submod_192r_3w", package = "dndscv")
# consensus_file = "/workspace/datasets/transfer/ferriol_deepcsa/test/expected_mutrate/consensus.exons_splice_sites.tsv"
# consensus_file = "/workspace/nobackup/bladder_ts/results/2024-11-03_deepCSA/createpanels/consensuspanels/consensus.exons_splice_sites.tsv"
# consensus_file = "/workspace/nobackup/bladder_ts/results/2024-11-03_deepCSA/createpanels/consensuspanels/consensus.exons_splice_sites.tsv"
# deepcsa_folder = "/workspace/datasets/transfer/ferriol_deepcsa/test/expected_mutrate/axel_test/"
# input_muts_file = "/workspace/nobackup/bladder_ts/results/2024-11-03_deepCSA/somaticmutations/all_samples.somatic.mutations.tsv"
# output_folder = "~/Documents/deepCSA_test"
# sample_depths_file = "/workspace/nobackup/bladder_ts/results/2024-11-03_deepCSA/annotatedepths/all_samples_indv.depths.tsv.gz"
# annotated_panelfile = "/workspace/nobackup/bladder_ts/results/2024-11-03_deepCSA/createpanels/panelannotation/captured_panel.tab.gz"
# input_muts_file = "/workspace/nobackup/bladder_ts/results/2024-11-03_deepCSA/somaticmutations/all_samples.somatic.mutations.tsv"
# output_folder = "~/Documents/deepCSA_test"
# sample_depths_file = "/workspace/nobackup/bladder_ts/results/2024-11-03_deepCSA/annotatedepths/all_samples_indv.depths.tsv.gz"
# annotated_panelfile = "/workspace/nobackup/bladder_ts/results/2024-11-03_deepCSA/createpanels/panelannotation/captured_panel.tab.gz"

# match dnds mutation types with the canonical mutation types
mut_types = rownames(substmodel)
# match dnds mutation types with the canonical mutation types
mut_types = rownames(substmodel)

consensus_file = args[1]
input_muts_file = args[2]
sample_depths_file = args[3]
annotated_panelfile = args[4]
output_folder = args[5]
input_file_intervals = args[6]

sample_sex_metadata_file = ifelse(length(args) >= 7, args[7], FALSE)
# sample_sex_metadata_file = "/data/bbg/projects/bladder_ts/data/complete_cohort/samples_metadata/complete_cohort_bladder.discarded_histofindings_lowmuts.clinical_variables_extended.no_transposed.tsv"
# TODO
# document this:
# this file must contain a row per sample and at least two columns called SAMPLE_ID and SEX, where male sex is encoded as "M"


# define functions:
inclusion_exclusion <- function(plist) {
  n <- length(plist)

  if (n > 2) {
    p1 <- inclusion_exclusion(plist[1:(n %/% 2)])
    p2 <- inclusion_exclusion(plist[(n %/% 2 + 1):n])
    return(inclusion_exclusion(c(p1, p2)))
  }

  if (n == 2) {
    return(plist[1] + plist[2] - plist[1] * plist[2])
  }

  if (n == 1) {
    return(plist[1])
  }
}


inclusion_exclusion <- function(plist) {
  n <- length(plist)

  if (n > 2) {
    p1 <- inclusion_exclusion(plist[1:(n %/% 2)])
    p2 <- inclusion_exclusion(plist[(n %/% 2 + 1):n])
    return(inclusion_exclusion(c(p1, p2)))
  }

  if (n == 2) {
    return(plist[1] + plist[2] - plist[1] * plist[2])
  }

  if (n == 1) {
    return(plist[1])
  }
}


get_dnds_mut_context = function(muts) {

  muts_strand = muts |>
    mutate(genestrand = ifelse(strand == 1, "+", "-"))


  muts_gr = GRanges(seqnames = muts_strand$chr, IRanges(muts_strand$pos, end = muts_strand$pos)) + 1
  strand(muts_gr) = Rle(muts_strand$genestrand)

  muts_strand = muts_strand |>
    mutate(ctx = getSeq(Hsapiens, muts_gr, as.character =  TRUE)) |>
    mutate(refstrand = case_when(genestrand == "-" & ref == "A" ~ "T",
                                 genestrand == "-" & ref == "C" ~ "G",
                                 genestrand == "-" & ref == "G" ~ "C" ,
                                 genestrand == "-" & ref == "T" ~ "A", .default = ref),
           altstrand = case_when(genestrand == "-" & alt == "A" ~ "T",
                                 genestrand == "-" & alt == "C" ~ "G",
                                 genestrand == "-" & alt == "G" ~ "C" ,
                                 genestrand == "-" & alt == "T" ~ "A", .default = alt)) |>
    mutate(mutation = paste0(ctx, ">", substr(ctx, 1,1), altstrand, substr(ctx, 3,3)))

  # collect all splicing mutations and gather them under the 'splicing' part
  muts_strand = muts_strand |>
    mutate(impact = ifelse(grepl("splice", IMPACT), "spl", IMPACT)) |>
    mutate(GENE = factor(GENE))

  return(muts_strand)
}

generate_tables = function(tab) {

  # initialize empty matrix
  mat = matrix(0, 192, 4, dimnames = list(mut_types, c("synonymous", "missense", "nonsense", "spl")))

  # initialize empty gene mutation list
  gene_tabs = list()
  for (i in levels(tab$GENE)) {
    gene_tabs[[i]] = as.data.frame(mat)
  }

  # strand: get gene strands
  gene_tables = tab |>
    group_by(GENE, impact, mutation) |>
    count() |>
    pivot_wider(names_from = impact, values_from = n, values_fill = 0) |>
    ungroup()

  gene_tables = split(gene_tables |> select(-GENE), gene_tables$GENE)
  gene_tables = lapply(gene_tables, column_to_rownames, "mutation")

  for (gene in names(gene_tables)) {
    rows = rownames(gene_tables[[gene]])
    columns = colnames(gene_tables[[gene]])
    gene_tabs[[gene]][rows, columns]  = gene_tables[[gene]]
  }

  gene_tabs = lapply(gene_tabs, as.matrix)
  return(gene_tabs)
}

generate_tables_depth = function(tab_depths, sample) {

  # initialize empty matrix
  mat = matrix(0, 192, 4, dimnames = list(mut_types,
   c("synonymous", "missense", "nonsense", "spl")))

  # initialize empty gene mutation list
  gene_tabs = list()
  for (i in levels(tab_depths$GENE)) {
    gene_tabs[[i]] =  as.data.frame(mat)
  }

  # correct for the depth of the sequenced sample. If no sample is selected, select all patient columns
  if (missing(sample)) {
    site_depths = rowSums(tab_depths |>
                            select(-c("chr", "pos", "ref", "alt",
                             "MUT_ID", "GENE", "IMPACT", "CONTEXT_MUT", "CONTEXT",
                              "strand", "genestrand", "ctx", "refstrand", "altstrand", "mutation", "impact")))
  }   else {
    site_depths = tab_depths |>
      pull(all_of(sample))
  }

  # strand: get gene strands
  gene_tables = tab_depths |>
    mutate(depth = site_depths) |>
    group_by(GENE, impact, mutation) |>
    summarize(n = sum(depth)) |>
    pivot_wider(names_from = impact, values_from = n, values_fill = 0) |>
    ungroup()

  gene_tables = split(gene_tables |> dplyr::select(-GENE), gene_tables$GENE)
  gene_tables = lapply(gene_tables, column_to_rownames, "mutation")

  for (gene in names(gene_tables)) {
    rows = rownames(gene_tables[[gene]])
    columns = colnames(gene_tables[[gene]])
    gene_tabs[[gene]][rows, columns]  = gene_tables[[gene]]
  }

  gene_tabs = lapply(gene_tabs, as.matrix)
  return(gene_tabs)
}

# Subfunction from dNdScv: fitting substitution model using custom matrices from the duplex consensus regions
fit_substmodel = function(N, L, substmodel) {

  l = c(L); n = c(N); r = c(substmodel)
  n = n[l!=0]; r = r[l!=0]; l = l[l!=0]

  params = unique(base::strsplit(x=paste(r,collapse="*"), split="\\*")[[1]])
  indmat = as.data.frame(array(0, dim=c(length(r),length(params))))
  colnames(indmat) = params
  for (j in 1:length(r)) {
    indmat[j, base::strsplit(r[j], split="\\*")[[1]]] = 1
  }

  model = glm(formula = n ~ offset(log(l)) + . -1, data=indmat, family=poisson(link=log))
  mle = exp(coefficients(model)) # Maximum-likelihood estimates for the rate params
  ci = exp(confint.default(model)) # Wald confidence intervals
  par = data.frame(name=gsub("\`","",rownames(ci)), mle=mle[rownames(ci)], cilow=ci[,1], cihigh=ci[,2])
  return(list(par=par, model=model))
}

dndsloc = function(genemuts, RefCDS, constrain_wnon_wspl = TRUE, onesided = FALSE) {
  # FIXME
  # this function fails when there is no mutation of a given gene-impact type
  # there are some divisions by 0
  message("[4] Running dNdSloc...")

  locll = function(nobs,nexp,x,indneut) {
    mrfold = max(1e-10, sum(nobs[indneut])/sum(nexp[indneut])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
    w = rep(1,4)
    w[-indneut] = nobs[-indneut]/nexp[-indneut]/mrfold
    w[nexp==0] = 0 # Suppressing cases where the expected rate is 0 (e.g. splice site mutations in genes with 1 exon)
    ll = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(w,dim=c(4,numrates))), log=T)) # loglik
    return(list(ll=ll,w=w))
  }

  lrt_onesided = function(ll0,llmis,lltrunc,llall,w) {
    # Under the 2 w model (wnon==wspl), one-sided tests can be performed using the log-likelihoods already calculated under different levels of constrain.
    # Positive selection: H0: wmis<=1, wtrunc<=1. H1: wmis>1 | wtrunc>1.
    # Negative selection: H0: wmis>=1, wtrunc>=1. H1: wmis<1 | wtrunc<1.
    if (w[1]>=1 & w[2]>=1) {
      ll0pos = ll0
      ll0neg = llall
    } else if (w[1]>=1 & w[2]<=1) {
      ll0pos = llmis
      ll0neg = lltrunc
    } else if (w[1]<=1 & w[2]>=1) {
      ll0pos = lltrunc
      ll0neg = llmis
    } else if (w[1]<=1 & w[2]<=1) {
      ll0pos = llall
      ll0neg = ll0
    }
    # One-sided LRTs: conservatively assuming 2 df even when under mixtures of positive and negative selection, the range of some parameters is constrained (1 df)
    p = 1-pchisq(2*(llall-c(ll0pos,ll0neg)),df=c(2,2))
    return(p=p)
  }

  selfun_loc = function(j) {
    y = as.numeric(genemuts[j,-1])
    nobs = y[1:4] # Number of observed mutations in the gene (Synonymous, Missense, Nonsense, Splice)
    nexp = y[5:8] # Number of expected mutations in the gene (Synonymous, Missense, Nonsense, Splice)
    x = RefCDS[[j]]

    # Alternative likelihood models constraining specific classes of mutations to be neutral
    ll0 = locll(nobs,nexp,x,1:4)$ll # Neutral model: wmis==1, wnon==1, wspl==1
    llmis = locll(nobs,nexp,x,c(1,2))$ll # Missense model: wmis==1, wnon free, wspl free
    lltrunc = locll(nobs,nexp,x,c(1,3,4))$ll # Truncation model: wmis free, wnon==1, wspl==1
    h = locll(nobs,nexp,x,1) # Fully unconstrained free selection model: wmis free, wnon free, wspl free
    llall_unc = h$ll
    wfree = h$w[-1]; wfree[wfree>1e4] = 1e4 # MLEs for the free w parameters (values higher than 1e4 will be set to 1e4)

    if (constrain_wnon_wspl == 0) { # Free selection model: free wmis, free wnon, free wspl (called llall_unc)

      # Models allowing free wnon or free wspl
      llnon = locll(nobs,nexp,x,c(1,3))$ll # Nonsense model: wmis free, wnon==1, wspl free
      llspl = locll(nobs,nexp,x,c(1,4))$ll # Splice model: wmis free, wnon free, wspl==1

      # LRTs: the free selection model is the alternative hypothesis, and the partially or fully constrained models are the nulls.
      p = 1-pchisq(2*(llall_unc-c(llmis,llnon,llspl,lltrunc,ll0)),df=c(1,1,1,2,3))

    } else { # Partially constrained free selection model: free wmis, free wnon==wspl (called llall)

      # Model for truncating substitutions forcing wnon==wspl
      mrfold = max(1e-10, sum(nobs[1])/sum(nexp[1])) # Correction factor of "t" based on the obs/exp of synonymous mutations in the gene
      wmisfree = nobs[2]/nexp[2]/mrfold; wmisfree[nexp[2]==0] = 0
      wtruncfree = sum(nobs[3:4])/sum(nexp[3:4])/mrfold; wtruncfree[sum(nexp[3:4])==0] = 0
      wfree = c(wmisfree,wtruncfree,wtruncfree)
      llall = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,wfree),dim=c(4,numrates))), log=T)) # loglik

      # LRTs: the free selection model is the alternative hypothesis, and the partially or fully constrained models are the nulls. The missense models (H0 and H1) is the same as for constrain_wnon_wspl == 0.
      p = 1-pchisq(2*c(llall_unc-llmis,llall-c(lltrunc,ll0)),df=c(1,1,2)) # Notice that for lltrunc and ll0, there is one fewer df when using wnon==wspl

      # Adding on-sided tests results if requested by the user
      if (onesided == T) {
        p_onesided = lrt_onesided(ll0,llmis,lltrunc,llall,wfree[1:2]) # ppos and pneg calculated by the lrt_onesided function
        p = c(p, p_onesided)
      }
    }
    return(c(wfree,p))
  }

  sel_loc = as.data.frame(t(sapply(1:nrow(genemuts), function(j) selfun_loc(j))))
  if (constrain_wnon_wspl == 0) {
    colnames(sel_loc) = c("wmis_loc","wnon_loc","wspl_loc","pmis_loc","pnon_loc","pspl_loc","ptrunc_loc","pall_loc")
    sel_loc$qmis_loc = p.adjust(sel_loc$pmis_loc, method="BH")
    sel_loc$qnon_loc = p.adjust(sel_loc$pnon_loc, method="BH")
    sel_loc$qspl_loc = p.adjust(sel_loc$pspl_loc, method="BH")
  } else {
    if (onesided == F) {
      colnames(sel_loc) = c("wmis_loc","wnon_loc","wspl_loc","pmis_loc","ptrunc_loc","pall_loc")
      sel_loc$qtrunc_loc = p.adjust(sel_loc$ptrunc_loc, method="BH")
      sel_loc$qall_loc = p.adjust(sel_loc$pall_loc, method="BH")
    } else {
      colnames(sel_loc) = c("wmis_loc","wnon_loc","wspl_loc","pmis_loc","ptrunc_loc","pall_loc", "ppos_loc", "pneg_loc")
      sel_loc$qtrunc_loc = p.adjust(sel_loc$ptrunc_loc, method="BH")
      sel_loc$qall_loc = p.adjust(sel_loc$pall_loc, method="BH")
      sel_loc$qpos_loc = p.adjust(sel_loc$ppos_loc, method="BH")
      sel_loc$qneg_loc = p.adjust(sel_loc$pneg_loc, method="BH")
    }
    sel_loc$qmis_loc = p.adjust(sel_loc$pmis_loc, method="BH")
  }
  sel_loc = cbind(genemuts[,1:5],sel_loc)
  sel_loc = sel_loc[order(sel_loc$pall_loc,sel_loc$pmis_loc,-sel_loc$wmis_loc),]
  return(sel_loc)
}

# estimate the mutation rates:
estimate_rates = function(mle_submodel, genemuts, RefCDS, relative_rates) {

  results_list = list()
  for (column in c("mle", "cilow", "cihigh")) {

    par = mle_submodel
    parmle = setNames(par[,column], par[,1])

    # how does the expected mutation load work:
    mutrates = sapply(substmodel[,1], \(x) prod(parmle[base::strsplit(x,split="\\*")[[1]]]))

    # Rle(mutrates, RefCDS[[1]]$L[,2])
    # inclusion_exclusion(Rle(mutrates, RefCDS[[1]]$L[,2]))

    genemuts_expected = t(sapply(RefCDS, \(x) colSums(x$L*mutrates)))
    colnames(genemuts_expected) = c("exp_syn", "exp_mis", "exp_non", "exp_spl")
    df = data.frame(gene_name = sapply(RefCDS, \(x) x$gene_name),
                    genemuts_expected) |>
      column_to_rownames("gene_name")
    results_list[[column]] = df
  }

  # if covariates are used, used the "exp_syn_cv" value to correct the mutation loads
  if (min(genemuts$n_syn) > 10) {
    message("local mutation rate is sufficiently high across all studied mutations, local synonymous rates will be used (dNdSloc)")
    ratio = genemuts$n_syn / genemuts$exp_syn
  } else if ("exp_syn_cv" %in% colnames(genemuts)) {
    message("covariate-based adjustment of the mutation rates will be used (dNdScv)")
    ratio = genemuts$exp_syn_cv / results_list$mle$exp_syn
  } else if (!missing(relative_rates)) {
    message("user supplied relative rates will be used")
    ratio = relative_rates
  }
  else {
    message("no covariates / not enough synonymous mutations -> switching to dNdSglobal mutation rate")
    ratio = 1 # do not convert the rates
  }

  results_list_cv = lapply(results_list, \(x) x*ratio)
  results_list_cv = lapply(results_list_cv, \(x) rownames_to_column(.data = x,var =  "gene_name"))

  mutation_estimates = rbindlist(results_list_cv, idcol = "type") |>
    pivot_longer(c(-type, -gene_name), names_to = "consequence") |>
    pivot_wider(names_from = "type", values_from = "value") |>
    dplyr::rename(mutrate = "mle")
  return(mutation_estimates)
}


# estimate the mutation rates with inclusion/exclusion principle:
estimate_rates_ie = function(mle_submodel, genemuts, RefCDS_1_genome, relative_rates) {

  # if covariates are used, used the "exp_syn_cv" value to correct the mutation loads
  if (min(genemuts$n_syn) > 10) {
    message("local mutation rate is sufficiently high across all studied mutations, local synonymous rates will be used (dNdSloc)")
    ratio = genemuts$n_syn / genemuts$exp_syn
  } else if ("exp_syn_cv" %in% colnames(genemuts)) {
    message("covariate-based adjustment of the mutation rates will be used (dNdScv)")
    ratio = genemuts$exp_syn_cv / results_list$mle$exp_syn
  } else if (!missing(relative_rates)) {
    message("user supplied relative rates will be used")
    ratio = relative_rates
  } else {
    message("no covariates / not enough synonymous mutations -> switching to dNdSglobal mutation rate")
    ratio = 1 # do not convert the rates
  }


  results_list = list()
  for (column in c("mle", "cilow", "cihigh")) {

    par = mle_submodel
    parmle = setNames(par[,column], par[,1])

    # how does the expected mutation load work:
    mutrates = sapply(substmodel[,1], \(x) prod(parmle[base::strsplit(x,split="\\*")[[1]]]))

    # make indexes for combinations. Trunc is nonsense + splicing mutations
    index_consequences = list(exp_syn = 1, exp_mis = 2, exp_non = 3, exp_spl = 4, exp_trunc = c(3, 4), exp_all = 1:4)
    df_exp = matrix(NA, nrow = length(RefCDS_1_genome), ncol = length(index_consequences))
    rownames(df_exp) = sapply(RefCDS, \(x) x$gene_name)
    colnames(df_exp) = names(index_consequences)
    df_exp_ie = df_exp

    # multiply the rates to the lenght of each of the coding sequence.
    for (consequence in names(index_consequences)) {
      index = index_consequences[[consequence]]
      consequence_rates = mapply(RefCDS_1_genome, ratio, FUN = \(x,y) {
             rates = Rle(mutrates, rowSums(x$L[,index, drop = FALSE])) |> as.numeric()
             if (length(rates) == 0 ) {rates = 0}
             rates = rates * y
             c(sum(rates), inclusion_exclusion(rates))})

      df_exp[,consequence] = consequence_rates[1,]
      df_exp_ie[,consequence] = consequence_rates[2,]
    }

    df_exp = as.data.frame(df_exp) |> rownames_to_column("gene_name")
    df_exp_ie = as.data.frame(df_exp_ie) |> rownames_to_column("gene_name")

    results_list[[column]] = rbindlist(list(sum = df_exp, inclusion_exclusion = df_exp_ie), idcol = "correction") |>
      pivot_longer(starts_with("exp"), names_to = "consequence", values_to = "mutrate")
  }

  mutation_estimates = rbindlist(results_list, idcol = "column") |>
    pivot_wider(values_from = "mutrate", names_from = "column")

  return(mutation_estimates)
}


# estimate the mutation rates with inclusion/exclusion principle:
estimate_rates_ie = function(mle_submodel, genemuts, RefCDS_1_genome, relative_rates) {

  # if covariates are used, used the "exp_syn_cv" value to correct the mutation loads
  if (min(genemuts$n_syn) > 10) {
    message("local mutation rate is sufficiently high across all studied mutations, local synonymous rates will be used (dNdSloc)")
    ratio = genemuts$n_syn / genemuts$exp_syn
  } else if ("exp_syn_cv" %in% colnames(genemuts)) {
    message("covariate-based adjustment of the mutation rates will be used (dNdScv)")
    ratio = genemuts$exp_syn_cv / results_list$mle$exp_syn
  } else if (!missing(relative_rates)) {
    message("user supplied relative rates will be used")
    ratio = relative_rates
  } else {
    message("no covariates / not enough synonymous mutations -> switching to dNdSglobal mutation rate")
    ratio = 1 # do not convert the rates
  }


  results_list = list()
  for (column in c("mle", "cilow", "cihigh")) {

    par = mle_submodel
    parmle = setNames(par[,column], par[,1])

    # how does the expected mutation load work:
    mutrates = sapply(substmodel[,1], \(x) prod(parmle[base::strsplit(x,split="\\*")[[1]]]))

    # make indexes for combinations. Trunc is nonsense + splicing mutations
    index_consequences = list(exp_syn = 1, exp_mis = 2, exp_non = 3, exp_spl = 4, exp_trunc = c(3, 4), exp_all = 1:4)
    df_exp = matrix(NA, nrow = length(RefCDS_1_genome), ncol = length(index_consequences))
    rownames(df_exp) = sapply(RefCDS, \(x) x$gene_name)
    colnames(df_exp) = names(index_consequences)
    df_exp_ie = df_exp

    # multiply the rates to the lenght of each of the coding sequence.
    for (consequence in names(index_consequences)) {
      index = index_consequences[[consequence]]
      consequence_rates = mapply(RefCDS_1_genome, ratio, FUN = \(x,y) {
             rates = Rle(mutrates, rowSums(x$L[,index, drop = FALSE])) |> as.numeric()
             if (length(rates) == 0 ) {rates = 0}
             rates = rates * y
             c(sum(rates), inclusion_exclusion(rates))})

      df_exp[,consequence] = consequence_rates[1,]
      df_exp_ie[,consequence] = consequence_rates[2,]
    }

    df_exp = as.data.frame(df_exp) |> rownames_to_column("gene_name")
    df_exp_ie = as.data.frame(df_exp_ie) |> rownames_to_column("gene_name")

    results_list[[column]] = rbindlist(list(sum = df_exp, inclusion_exclusion = df_exp_ie), idcol = "correction") |>
      pivot_longer(starts_with("exp"), names_to = "consequence", values_to = "mutrate")
  }

  mutation_estimates = rbindlist(results_list, idcol = "column") |>
    pivot_wider(values_from = "mutrate", names_from = "column")

  return(mutation_estimates)
}

##### LOAD DATA #######
# Read in consensus file and mutations file:
consensus = fread(consensus_file) |>
  dplyr::rename(chr = CHROM, pos = POS, ref = REF, alt = ALT) |>
  mutate(GENE = as.factor(GENE))

input_muts = fread(input_muts_file) |>
  select(CHROM, POS, REF, ALT, SYMBOL, STRAND, SAMPLE_ID) |>
  dplyr::rename(sampleID = SAMPLE_ID, chr = CHROM, pos = POS, ref = REF, alt = ALT, strand = STRAND, GENE = SYMBOL) |>
  select(sampleID, chr, pos, ref, alt, strand, GENE)

# read in the depth file table for all samples:
sample_depths = fread(sample_depths_file)
sample_depths =  sample_depths |>
  dplyr::select(-CONTEXT) |>
  dplyr::rename(chr = CHROM, pos = POS)


### ANALYSIS ####
consensus_depths = left_join(consensus, sample_depths, by = c("chr", "pos"))

# alternative gene strands approach:
command = sprintf("zcat % s | cut -f16,18 | grep -v '##' | grep -v 'SYMBOL' | sort | uniq", annotated_panelfile)
gene_strands = read.table(text = system(command, intern = TRUE), col.names = c("strand", "GENE"))

consensus_depths = left_join(consensus_depths,gene_strands) |>
  filter(!is.na(strand))
consensus_depths[is.na(consensus_depths)] = 0

# join mutations with consensus information
muts = inner_join(input_muts, consensus)
Nlist = muts |> get_dnds_mut_context() |> generate_tables()
Llist = consensus_depths |>  get_dnds_mut_context() |> generate_tables_depth()

# bind and sum the matrix to perform Poisson regression
Nall = abind::abind(Nlist, along = 3)
Lall = abind::abind(Llist, along = 3)
N = apply(Nall, c(1,2), sum)
L = apply(Lall, c(1,2), sum)

# Fitting all mutation rates and the 3 global selection parameters
poissout = fit_substmodel(N, L, substmodel) # Original substitution model
par = poissout$par
poissmodel = poissout$model
parmle =  setNames(par[,2], par[,1])
mle_submodel = par
rownames(mle_submodel) = NULL

# Fitting models with 1 and 2 global selection parameters
s1 = gsub("wmis","wall",gsub("wnon","wall",gsub("wspl","wall",substmodel)))
par1 = fit_substmodel(N, L, s1)$par # Substitution model with 1 selection parameter
s2 = gsub("wnon","wtru",gsub("wspl","wtru",substmodel))
par2 = fit_substmodel(N, L, s2)$par # Substitution model with 2 selection parameter
globaldnds = rbind(par, par1, par2)[c("wmis","wnon","wspl","wtru","wall"),]
sel_loc = sel_cv = NULL

## 4. dNdSloc: variable rate dN/dS model (gene mutation rate inferred from synonymous subs in the gene only)
genemuts = data.frame(gene_name = names(Nlist), n_syn=NA, n_mis=NA, n_non=NA, n_spl=NA, exp_syn=NA, exp_mis=NA, exp_non=NA, exp_spl=NA, stringsAsFactors=F)
genemuts[,2:5] = t(colSums(Nall))
mutrates = sapply(substmodel[,1], function(x) prod(parmle[base::strsplit(x,split="\\*")[[1]]])) # Expected rate per available site
genemuts[,6:9] = t(sapply(Llist, function(x) colSums(x*mutrates)))
numrates = length(mutrates)

# plot the number of expected mutations in the genes according to the model:
expected_genemuts = genemuts |>
  mutate(ratio = n_syn/exp_syn) |>
  mutate(across(starts_with("exp"), ~ .*ratio)) |>
  select(-ratio)

output_list$expected_genemuts = expected_genemuts

# # plot the ratio's for the mutation rate across genes:
plot_list$genemuts_ratio = genemuts |>
  mutate(ratio = n_syn / exp_syn)  |>
  ggplot(aes(x = gene_name, y = ratio)) +
  geom_col() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  labs(y = "Observed / Expected synonymous mutations", x = NULL,
       title = "Observed synonymous vs Expectd synonymous", subtitle = "Indicates the baseline mutation load of a gene") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# make sure the RefCDS in this case is the list of duplex-sequenced genes
RefCDS = list()
for (i in 1:length(Nlist)) {
  name = names(Nlist)[i]
  RefCDS[[i]] = list(L = Llist[[i]], N = Nlist[[i]], gene_name = name  )
}


# run dNdSloc across genes
sel_loc_trunc = dndsloc(genemuts, RefCDS, constrain_wnon_wspl = TRUE, onesided = FALSE)
sel_loc_all = dndsloc(genemuts, RefCDS, constrain_wnon_wspl = FALSE, onesided = FALSE)

rownames(sel_loc_trunc) = sel_loc_trunc$gene_name

dndstrunc_heatmap = pheatmap::pheatmap(sel_loc_trunc[,10:14], angle_col = 45, display_numbers = TRUE, main = "P-values dNdSloc \ngene regions matched to consensus area duplex capture")
png(paste0(output_folder, "/dndstrunc_heatmap.png"), width = 1800, height = 1800, res = 250)
dndstrunc_heatmap
dev.off()

rownames(sel_loc_all) = sel_loc_all$gene_name
dndsall_heatmap = pheatmap::pheatmap(sel_loc_all[,10:14], angle_col = 45, display_numbers = TRUE, main = "P-values dNdSloc \ngene regions matched to consensus area duplex capture")
png(paste0(output_folder, "/dndsall_heatmap.png"), width = 1800, height = 1800, res = 250)
dndsall_heatmap
dev.off()

# plot of the dNdS values for the truncating mutations
sel_loc_trunc |>
  dplyr::select(gene_name, wmis_loc, wnon_loc) |>
  dplyr::rename(wtrunc_loc = wnon_loc) |>
  pivot_longer(-gene_name) |>
  ggplot(aes(x = gene_name, y = value, fill = name)) + geom_col(position = "dodge") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(fill = NULL, y = "dNdS w value", x = NULL, title = "dNdSloc duplex data using deepCSA consensus regions") +
  theme_bw() +
  ggsci::scale_fill_nejm()
ggsave(paste0(output_folder, "/dnds_bar_trunc.png"), width = 11, height = 4)

# plot of the non-truncating mutations
dNdS_bar_all = sel_loc_all |>
  dplyr::select(gene_name, wmis_loc, wnon_loc, wspl_loc) |>
  pivot_longer(-gene_name) |>
  ggplot(aes(x = gene_name, y = value, fill = name)) + geom_col(position = "dodge") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(fill = NULL, y = "dNdS w value", x = NULL, title = "dNdSloc duplex data using deepCSA consensus regions") +
  theme_bw() +
  ggsci::scale_fill_nejm()
ggsave(paste0(output_folder, "/dnds_bar_all.png"), dNdS_bar_all, width = 11, height = 4)

### Model the mutation rate for 1 single genome:
# make sure the RefCDS does contain a genome length of '1' genome, instead of the depth correction:
Llist = consensus_depths |>  get_dnds_mut_context() |> generate_tables()
Lall = abind::abind(Llist, along = 3)
L = apply(Lall, c(1,2), sum)

# make sure the RefCDS in this case is the list of duplex-sequenced genes
RefCDS_1_genome = list()
for (i in 1:length(Llist)) {
  RefCDS_1_genome[[i]] = list(L = Llist[[i]], N = Nlist[[i]], gene_name = names(Llist)[i] )
}

# calculate the total expected number of mutations
mutation_rates_avg = estimate_rates_ie(mle_submodel, genemuts = genemuts, RefCDS_1_genome = RefCDS_1_genome)

# intergrate with chromosome number
mutation_rates_avg
gene_chrs = consensus_depths |>
  select(chr, GENE) |> distinct() |>
  dplyr::rename(gene_name = GENE)
per_cell = left_join(gene_chrs, mutation_rates_avg) |>
  mutate(across(c(mle, cilow, cihigh), \(x) x*2))

# plot the difference across methods to determine the overall mutated epithelium.
per_cell |>
  filter(consequence == "exp_all") |>
  ggplot(aes(x = gene_name, y = mle, fill = correction)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = cilow, ymax = cihigh), position = position_dodge(width = 0.9), width = 0.2) +
  theme_classic() +
  labs(y = "fraction of mutated cells") +
  scale_y_continuous(expand=expansion(mult=c(0,0.1)))

ggsave(paste0(output_folder, "/compare_effect_inclusion_exclusion_average_cohort.png"),  width = 10, height = 6)
write.table(mutation_rates_avg, paste0(output_folder, "/mutation_rates_avg.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Estimation of the mutated fraction in all the mutated samples individually
# multiply all the samples to the mutation rate, mutational profile and the depth of each of the samples.
muts_context = muts |>  get_dnds_mut_context()
total_syn = muts_context |>
  filter(impact == "synonymous") |>
  group_by(mutation) |> count()

sample_syn = muts_context |>
  filter(impact == "synonymous") |>
  group_by(mutation, sampleID) |>
  count() |>
  pivot_wider(names_from = sampleID, values_from = n, values_fill = 0) |>
  ungroup()
samples_syn_rel = sample_syn |>
  mutate(across(-mutation, ~ . /
                  total_syn$n))

# Using all mutations in the model, we are better at handling sparse data (zero's)
consensus_depths_ctx = consensus_depths |>  get_dnds_mut_context()

rate_list = list()
globaldnds_list = list()

# calculate the 'relative mutation rate across the cohort, which can then be applied to each sample individually
ratios_all_samples = genemuts$n_syn / genemuts$exp_syn
names(ratios_all_samples) = genemuts$gene_name
names(ratios_all_samples) = genemuts$gene_name

for (sample_name in unique(muts$sampleID)) {
  # check mutation rates using dnds:
  print(sample_name)
  Nlist = muts_context |>
    filter(sampleID == sample_name) |>
    generate_tables()

  Llist = consensus_depths_ctx |> generate_tables_depth(sample = sample_name)
  Lall = abind::abind(Llist, along = 3)
  Nall = abind::abind(Nlist, along = 3)
  L = apply(Lall, c(1,2), sum)
  N = apply(Nall, c(1,2), sum)

  # Fitting all mutation rates and the 3 global selection parameters
  poissout = fit_substmodel(N, L, substmodel) # Original substitution model
  par = poissout$par
  poissmodel = poissout$model
  parmle =  setNames(par[,2], par[,1])
  mle_submodel = par
  rownames(mle_submodel) = NULL

  # Fitting models with 1 and 2 global selection parameters
  s1 = gsub("wmis","wall",gsub("wnon","wall",gsub("wspl","wall",substmodel)))
  par1 = fit_substmodel(N, L, s1)$par # Substitution model with 1 selection parameter
  s2 = gsub("wnon","wtru",gsub("wspl","wtru",substmodel))
  par2 = fit_substmodel(N, L, s2)$par # Substitution model with 2 selection parameter
  globaldnds = rbind(par, par1, par2)[c("wmis","wnon","wspl","wtru","wall"),]
  sel_loc = sel_cv = NULL # selection values for the single samples

  ## 4. dNdSloc: variable rate dN/dS model (gene mutation rate inferred from synonymous subs in the gene only)
  genemuts_sample = data.frame(gene_name = names(Llist), n_syn=NA, n_mis=NA, n_non=NA, n_spl=NA, exp_syn=NA, exp_mis=NA, exp_non=NA, exp_spl=NA, stringsAsFactors=F)
  genemuts_sample[,2:5] = t(colSums(Nall))
  mutrates = sapply(substmodel[,1], function(x) prod(parmle[base::strsplit(x,split="\\*")[[1]]])) # Expected rate per available site
  genemuts_sample[,6:9] = t(sapply(Llist, function(x) colSums(x*mutrates)))

  # estimate the mutation rates
  mutation_rates_1_genome = estimate_rates_ie(mle_submodel,genemuts = genemuts_sample, RefCDS_1_genome = RefCDS_1_genome, relative_rates = ratios_all_samples)
  mutation_rate_depths = mutation_rates_1_genome |>
    mutate(`expected mutations` = case_match(consequence,
                                             "exp_syn" ~ "synonymous",
                                             "exp_non" ~ "nonsense",
                                             "exp_mis" ~ "missense",
                                             "exp_spl" ~ "splicing",
                                             "exp_trunc" ~ "truncating"))

  rate_list[[sample_name]] = mutation_rate_depths
  globaldnds_list[[sample_name]] = globaldnds
}

fraction_genemut_genome = rbindlist(rate_list, idcol =  "sample")
output_list$fraction_genemut_genome = fraction_genemut_genome |>
    select(-sample, -`expected mutations`)

globaldnds_all = rbindlist(globaldnds_list, idcol =  "sample")
#plot the global dnds by sample.
dnds_by_sample = ggplot(globaldnds_all, aes(x = reorder(sample, mle), y = mle)) +
    facet_grid(name ~ . , scales = "free_y") +
    geom_col() +
    geom_errorbar(aes(ymin = cilow, ymax = cihigh), width = 0.5) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_hline(yintercept = 1, linetype = "dotted") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
    scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
    labs(x = NULL)
ggsave(paste0(output_folder, "/dnds_by_sample.png"), dnds_by_sample, width = 5 + (length(globaldnds_all)/3), height = 5)



# load metadata to know sex of patients
if (sample_sex_metadata_file != FALSE){
    print("Using sample metadata information")

    metadata_file = sample_sex_metadata_file
    metadata = fread(metadata_file) |>
        select(SAMPLE_ID, SEX) |>
        dplyr::rename(sample = SAMPLE_ID, sex = SEX) |>
        distinct()

    fraction_genemut_cell = left_join(fraction_genemut_genome, metadata)  |>
        left_join(gene_chrs) |>
        mutate(across(c(mle, cihigh, cilow), ~ ifelse(chr == "chrX" & sex == "M",  ., . * 2))) |>
        select(-`expected mutations`, -sex, -chr)

    output_list$fraction_genemut_cell = fraction_genemut_cell

    # testing plot - check the effect of inclusion-exclusion principle
    sample_cell_ie = fraction_genemut_cell |>
        filter(correction == "inclusion_exclusion" & consequence == "exp_all") |>
        group_by(sample) |>
        summarize(inclusion_exclusion = inclusion_exclusion(mle),
                    simple_sum = sum(mle)) |>
        arrange(desc(inclusion_exclusion)) |>
        mutate(sample = factor(sample, levels = sample)) |>
        pivot_longer(c(inclusion_exclusion, simple_sum))

    ggplot(sample_cell_ie, aes(x = sample, y = value, fill = name)) +
        geom_col(position = "dodge") +
        theme(axis.text.x = element_text(angle = 90)) +
        labs(x = NULL, y = "fraction of mutated genomes",
            title = "fraction of mutated genomes") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_y_continuous(expand=expansion(mult=c(0,0.1)))
    ggsave(paste0(output_folder, "/total_samples_difference_inclusion_exclusion.png"), width = 12, height = 6, dpi = 300)

    # values ready for publication:
    output_list$percentage_mutated_epithelium = sample_cell_ie |>
        pivot_wider() |>
        select(-simple_sum)

}

# save output tables
if (exists("output_folder")) {
    for (i in names(output_list)) {
        print(i)
        fwrite(output_list[[i]], paste0(output_folder, "/", i, ".tsv"), sep = "\t")
        }
    }



# save RefCDS object
saveRDS(RefCDS_1_genome, paste0(output_folder, "/RefCDS.rds"))

# Create a new object named RefCDS
RefCDS <- RefCDS_1_genome



#### Create gr_genes object needed for running dNdScv original package
df_intervals <- fread(input_file_intervals)

# Create the GRanges object
# Ensure column names are mapped correctly: CHROMOSOME -> seqnames, START/END -> ranges, ELEMENT -> gene names
gr_genes <- GRanges(
  seqnames = df_intervals$CHROMOSOME,
  ranges = IRanges(start = df_intervals$START, end = df_intervals$END),
  strand = "*"
)
# Assign gene names from the ELEMENT column
mcols(gr_genes)$names <- df_intervals$ELEMENT

# Save RefCDS and gr_genes together in the .rda file
save(RefCDS, gr_genes, file = file.path(output_folder, "RefCDS.rda"))
