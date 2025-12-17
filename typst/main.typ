#import "@preview/clean-math-paper:0.2.5": *

#let date = datetime.today().display("[month repr:long] [day], [year]")

// Modify some arguments, which can be overwritten in the template call
#page-args.insert("numbering", "1/1")
#text-args-title.insert("size", 2em)
#text-args-title.insert("fill", black)
#text-args-authors.insert("size", 12pt)

#show: template.with(
title: "Optimal selection of antisense oligonucleotides for haplotype-selective knockdown in diverse populations",
authors: (
(name: "Barney Hill", affiliation-id: "1,2"),
(name: "Duncan S. Palmer", affiliation-id: "2"),
(name: "Stephan J. Sanders", affiliation-id: "1,3,4,5"),
),
affiliations: (
(id: 1, name: "Department of Paediatrics, University of Oxford, OX3 7TY Oxford, United Kingdom"),
(id: 2, name: "Big Data Institute, Li Ka Shing Centre for Health Information and Discovery, University of Oxford, Oxford, United Kingdom"),
(id: 3, name: "Institute of Developmental and Regenerative Medicine, IMS-Tetsuya Nakamura Building, Old Road Campus, OX3 7TY Oxford, United Kingdom"),
(id: 4, name: "New York Genome Center, New York, NY 10013, USA"),
(id: 5, name: "Department of Psychiatry and Behavioral Sciences, UCSF Weill Institute for Neurosciences, University of California San Francisco, San Francisco, CA 94178, USA"),
(id: "*", name: "Correspondence to barney.hill@merton.ox.ac.uk")
),
date: date,
link-color: rgb("#008002"),
abstract: [Antisense oligonucleotides (ASOs) offer a promising therapeutic approach for gain-of-function _SCN2A_ mutations causing developmental and epileptic encephalopathy. Because _SCN2A_ is haploinsufficient, selective knockdown of only the disorder-bearing haplotype is preferable. ASOs can achieve this selectivity by targeting common variants that are heterozygous in a given patient and in cis with the pathogenic variant. However, when the variant-bearing haplotype is unknown a priori, selecting optimal ASO targets to maximise population coverage becomes a combinatorial problem. We formalise this optimisation for both single-variant selection and multi-variant portfolios to allow for redundancy during preclinical development. Applying our framework to 4,094 globally diverse phased genomes, we show that five optimally selected indels in _SCN2A_ achieve approximately 80% expected population coverage, with equitable performance across ancestral groups. This approach provides a principled strategy for designing allele-selective ASO programmes that can serve broad patient populations.],
)

= Introduction

== Clinical Background

- _SCN2A_-associated developmental and epileptic encephalopathy (DEE) is a rare monogenic disorder characterised by early-onset seizures and profound developmental impairment
- Gain-of-function (GoF) mutations cause neuronal hyperexcitability, making knockdown of the mutant allele a rational therapeutic strategy
- _SCN2A_ is haploinsufficient: loss of the wild-type allele causes a distinct loss-of-function phenotype, making preservation of wild-type expression essential

== ASO Strategies for Haploinsufficient Genes

Antisense oligonucleotides offer two broad approaches for GoF mutations in haploinsufficient genes:

*Non-allele-selective ASOs* use a single design to treat all patients regardless of their specific mutation. Previously there have been reported candidates with promising preclinical efficacy @li_antisense_2021 and recently Elsunersen, targeting the _SCN2A_ transcript non-selectively, achieved >60% seizure reduction in a GoF patient @wagner_antisense_2025. However, this approach risks reducing wild-type expression below the haploinsufficiency threshold, posing toxicity concerns.

*Allele-selective ASOs* preserve wild-type expression by targeting only the mutation-bearing haplotype. The challenge lies in selecting appropriate targets: each patient's unique de novo mutation would conventionally require bespoke ASO design, limiting scalability.

== Targeting Common Variants in Cis

A scalable allele-selective strategy targets common variants that are heterozygous in the patient and in cis with the pathogenic mutation. nLorem has recently demonstrated this approach for _SCN2A_, disclosing two gapmers targeting benign SNPs in cis with GoF mutations @kim-mcmanus_individualized_2025 @kingsmore_diplotype-based_2025. After screening over 500 candidate ASOs, they have administered (intrathecal) individualised treatments to two patients based on their compatible SNP genotypes.

#figure(
  kind: table,
  caption: [
    nLorem allele-selective ASOs for _SCN2A_.
  ],
)[
  #table(
    columns: 4,
    align: (left, left, left, left),
    stroke: 0.5pt,
    inset: 8pt,

    table.header(
      [*Patient*],
      [*Target SNP*],
      [*Allele freq*],
      [*Sequence*],
    ),

    [Patient 1],
    [rs190030016],
    [0.0002],
    [#text(fill: rgb("#1E90FF"))[TGoCoCoAo]ACAATGTAC#text(fill: rgb("#1E90FF"))[AoGoGGT]],

    [Patient 2],
    [rs72874313],
    [0.1789],
    [#text(fill: rgb("#1E90FF"))[GGoToCoAo]ATTGAAAGATAo#text(fill: rgb("#1E90FF"))[ToCCC]],
  )
]

#v(0.8em)
#text(size: 9pt)[
  *Key for ASO chemical modifications:*  
  black = unmodified deoxyribose;  
  #text(fill: rgb("#1E90FF"))[blue] = 2′-methoxyethyl (MOE).  
  Unmarked backbone linkages = phosphorothioate (PS);  
  linkages marked with o = normal phosphodiester (PO).  
  All cytosines have 5-methyl modifications (5-methylcytosine).
]

A second paper presents the use of Illumina's synthetic long-read sequencing platform for resolving patient haplotypes @cheng_constellation_2025. A screen of infants with seizure-related disorders found 16% were compatible with the two existing ASOs.

Additionally nLorem have tested this approach in the case of KIF1A @ziegler_antisense_2024.

== This Work

Here, we formalise this SNP optimisation problem: given candidate indels heterozygous across a population, select the $k$ targets that maximise expected coverage under uncertainty about which haplotype carries each patient's de novo mutation. We extend this framework to portfolios of redundant ASO candidates targeting variants in high linkage disequilibrium, mitigating the risk of individual candidates failing during preclinical development.

Applying our framework to 4,094 globally diverse phased genomes, we show that five optimally selected indels in _SCN2A_ achieve approximately 80% expected population coverage, with relatively equitable performance across ancestral groups. The FDA and MHRA is considering new guidance on basket and umbrella trials for rare disease therapies @prasad_fdas_nodate @noauthor_rare_nodate, opening the possibility of parallel trials for multiple ASOs matched to patient genetics.

// underselling the work - it is not just assessing the combinatorics, but creating a method to prioritize the SNPs for the nth gene

= Model

== Problem formulation

We consider a population of $M$ individuals with phased genotypes across $V$ candidate alleles. For each patient, the pathogenic mutation resides on one parental haplotype, but which one is unknown a priori. An ASO targeting allele $v$ achieves selective silencing only if the patient carries $v$ on the pathogenic haplotype but not on the wild-type haplotype. Given this constraint and the uncertainty about mutation location, we seek to select $k$ target alleles (or portfolios of alleles in high LD) to maximise expected population coverage.

== Allele-selective targeting requires heterozygosity

An ASO targeting allele $v$ can selectively silence one haplotype only if the patient carries $v$ on that haplotype but not on the other. Homozygosity for the target allele would result in knockdown of both haplotypes, which is undesirable for a haploinsufficient gene.

We encode phased genotypes as
$ H_(i v h) = cases(1 &"if individual" i "carries allele" v "on haplotype" h, 0 &"otherwise") $
for individuals $i in {1,...,M}$, alleles $v in {1,...,V}$, and haplotypes $h in {1,2}$.

We further define the heterozygosity indicator
$ tilde(H)_(i v h) = H_(i v h) (1 - H_(i v overline(h))) $
where $overline(h)$ denotes the opposite haplotype. This quantity equals 1 only when individual $i$ carries allele $v$ on haplotype $h$ but not on haplotype $overline(h)$.


== Modelling de novo mutations at the population level

Our goal is to select ASO targets that will serve a broad patient population by target a common variant rather then the specific pathogenic variant. Since germline de novo mutations arise randomly on either parental haplotype with equal probability, we model the disease-bearing haplotype as $h in {1, 2}$ with probability $1/2$ each.

== The single-ASO case

Consider first the simplest scenario in which we develop one ASO targeting allele $v$ and ask whether it can treat individual $i$. If the mutation were known to reside on haplotype $h$, the patient would be treatable if, and only if, $tilde(H)_(i v h) = 1$ that is, allele $v$ is present on the disease haplotype and absent from the wild-type haplotype.

Let $S = {v}$ denote the selection of a single target allele. The expected coverage for individual $i$, marginalising over the unknown mutation-bearing haplotype, is

$ E["cov"_i | S] = 1/2 tilde(H)_(i v 1) + 1/2 tilde(H)_(i v 2). $

This corresponds to the case $N = K = 1$ in the general formulation below.

== Modelling preclinical attrition

ASO development has substantial attrition, with many candidates failing toxicity, specificity, or efficacy thresholds before reaching patients. To mitigate this risk, a programme may advance $N$ candidates into preclinical development, anticipating that only $K <= N$ will survive to clinical use.

Let $S$ denote a selection of $N$ candidate alleles. For individual $i$ on haplotype $h$, define

$ n_(i h) = sum_(v in S) tilde(H)_(i v h) $

as the count of selected candidates for which individual $i$ is heterozygous on haplotype $h$. The expected coverage—averaging over all $binom(N, K)$ possible subsets of $K$ survivors—is

$ E["cov"_(i h) | S] = 1 - binom(N - n_(i h), K) / binom(N, K). $

Marginalising over the unknown mutation-bearing haplotype:

$ E["cov"_i | S] = 1/2 sum_(h in {1,2}) [ 1 - binom(N - n_(i h), K) / binom(N, K) ]. $

We seek the selection $S^*$ maximising expected population coverage:

$ S^* = arg max_(S subset.eq {1,...,V}, |S| = N) sum_(i=1)^M E["cov"_i | S]. $

Exhaustive enumeration over all $binom(V, N)$ combinations is intractable for realistic values of $V$ and $N$. However, the problem can be solved exactly using integer linear programming (ILP). A function $f$ is concave if each unit increase in its argument yields a smaller gain than the previous: $f(n+1) - f(n)$ is decreasing. Concave maximisation problems over linear constraints admit exact solutions via ILP. Our coverage function $f(n) = 1 - binom(N-n, K) / binom(N, K)$ is concave: the marginal coverage gain from the $(n+1)$-th helpful candidate is

$ f(n+1) - f(n) = binom(N - n - 1, K - 1) / binom(N, K) $

which decreases in $n$. The case $N = K$ (all selected candidates survive) corresponds to the single-variant selection problem and is handled identically.

= Results (_SCN2A_)

To model diverse populations we use 4,094 phased whole genomes from the Human Genome Diversity Project and the 1000 Genomes Project @koenig_harmonized_2024. For transcript definitions we use Ensembl 110 @dyer_ensembl_2025. For a given candidate allele $v$, we make the following filters:

1. $0.05 < "freq"_v < 0.95$.
2. There is no variant $u$ with $0.01 < "freq"_u < 0.99$ in the interval $["pos"_v - 20, "pos"_v + 20]$.
3. Let $u$ be the most common allele at the site where $u!=v$. Then $|"len"_v - "len"_u| >= 2$ bp.
4. There is at least one candidate design fully overlapping the target site with OligoAI score $>= 50$ (predicted 50% knockdown) and with no off-targets at $<= 1$ mismatch in hg38. * - NOT CURRENTLY IMPLEMENTED - TODO* 

Applying the implemented filters sequentially reduces the _SCN2A_ locus from 9,046 alleles to 50 candidates.

#let filter_counts = (
  ([Stage], [Filter], [Alleles remaining]),
  ([0], [All alleles within the _SCN2A_ locus window], [9,046]),
  ([1], [Minor-allele frequency 0.05–0.95], [896]),
  ([2], [No other common variant within $plus.minus 20$ bp], [746]),
  ([3], [$|"len"_v - "len"_u| >= 2$ bp vs. most common allele], [50]),
)

#table(
  columns: (auto, 1fr, auto),
  align: horizon + left,
  inset: 6pt,
  stroke: (x, y) => (
    bottom: if y == 0 { 1pt } else { 0.5pt + rgb("#dddddd") }
  ),
  ..filter_counts.flatten()
)

== Single-variant selection (P=1)

We first consider the case where each selected target corresponds to a single ASO (portfolio size $P = 1$). For each value of $k in {1, ..., 7}$, we perform exhaustive enumeration over all $binom(V, k)$ combinations of candidate variants to identify the subset $S^*$ maximising expected population coverage.

@fig-coverage-overall shows the maximum expected coverage as a function of the number of selected variants $k$. With a single optimally chosen indel ($k = 1$), we achieve approximately 25% expected coverage across the global population. Coverage increases monotonically with $k$, reaching approximately 80% at $k = 5$. The marginal gain diminishes as $k$ increases, reflecting the decreasing pool of uncovered individuals who are heterozygous at remaining candidate sites.

#figure(
  image("plots/aso_coverage_SCN2A_P=1_overall.png", width: 80%),
  caption: [Expected population coverage as a function of the number of selected ASO targets ($k$) with portfolio size $P = 1$. Coverage represents the probability that a randomly selected patient can be treated, averaged over uncertainty in which haplotype carries the pathogenic mutation.]
) <fig-coverage-overall>

== Ancestry-stratified coverage

Given the global diversity of our reference panel, we examined whether optimal variant selection provides equitable coverage across ancestral superpopulations. @fig-coverage-ancestry shows coverage stratified by the superpopulations: African (AFR), American (AMR), Central/South Asian ancestry (CSA), East Asian (EAS), European (EUR), Middle East (MID), and Oceanic (OCE).

Coverage is relatively uniform across populations except for individuals with Oceanic ancestry. Low coverage in Oceanic ancestry may be due to poor sample sizes (see N=28).

#figure(
  image("plots/aso_coverage_SCN2A_P=1_by_population.png", width: 80%),
  caption: [Expected coverage stratified by ancestral superpopulation for portfolio size $P = 1$. AFR: African; AMR: American; EAS: East Asian; EUR: European; SAS: South Asian. Disparities reflect population-specific allele frequency distributions.]
) <fig-coverage-ancestry>

To visualise the above results we can also display the selected alleles for $K=5, P=1$ and their coverage across a sample of individuals. Note how some selected alleles are reference (165340317_G_G) while others are alternative (165375368_A_AGT). Below we are displaying 

#let data = csv("plots/aso_haplotypes_k=5.csv")
#let header = data.first()
#let rows = data.slice(1, 15)
#let num_alleles=5


#table(
  columns: (auto, auto, ..((1fr,) * num_alleles), auto),
  align: center + horizon,
  stroke: (x, y) => (
    bottom: if y == 0 { 0.5pt } else if y == 1 { 1pt } else { 0.5pt + gray }
  ),
  
  // Top header row
  table.cell(rowspan: 2, align: horizon, strong[sample]), 
  table.cell(rowspan: 2, align: horizon, strong[superpop]), 
  table.cell(colspan: num_alleles, strong[Selected alleles (k=5)]),
  table.cell(rowspan: 2, align: horizon, strong[Coverage]),
  
  // Rotated allele headers only
  ..header.slice(2, -1).map(h => rotate(-90deg, reflow: true, text(size: 8pt, h))),
  
  // Data rows
  ..rows.flatten(),
  
  // Ellipsis row
  ..("...",) * header.len()
)

== Effect of portfolio size on coverage

To assess the impact of developing multiple ASO candidates per target site, we repeated the optimisation for portfolio sizes $P in {1,2,3,4,5}$. Recall that a portfolio of size $P$ contains the anchor variant and its $P - 1$ nearest neighbours by linkage disequilibrium, with one ASO from each portfolio assumed to survive preclinical development.

@fig-coverage-portfolio shows expected coverage across all populations as a function of $k$ for each portfolio size. Increasing $P$ provides substantial gains, particularly at lower values of $k$. For example, at $k = 3$, coverage increases from approximately 50% with $P = 1$ to approximately 65% with $P = 15$. This improvement arises because larger portfolios increase the probability that at least one ASO targets an allele present on the disease haplotype for any given patient.

#figure(
  image("plots/aso_coverage_SCN2A_portfolio_comparison.png", width: 80%),
  caption: [Expected population coverage for varying portfolio sizes $P in {1,2,3,4,5}$. Larger portfolios mitigate against preclinical failure by developing multiple redundant ASO candidates targeting variants in high linkage disequilibrium.]
) <fig-coverage-portfolio>

= Code Availability

All code for the optimisation framework and analysis is available at https://github.com/barneyhill/allele-specific.

#bibliography("zotero.bib")

= Supplementary Material

== Portfolio model for correlated candidates

An alternative approach to modelling preclinical attrition groups candidates into portfolios based on linkage disequilibrium. For each selected anchor allele $v$, we define a portfolio $cal(P)(v)$ containing $v$ and its $P - 1$ nearest neighbours by LD ($R^2$), with exactly one ASO per portfolio assumed to survive development.

Let $n_(i h p) = sum_(u in cal(P)(v_p)) tilde(H)_(i u h)$ count heterozygous alleles in portfolio $p$ for individual $i$ on haplotype $h$. The expected coverage for a selection of $k$ anchor alleles $S = {v_1, ..., v_k}$ is

$ E["cov"_i | S] = 1/2 sum_(h in {1,2}) [ 1 - product_(p=1)^k (1 - n_(i h p) / P) ]. $

The product over portfolios introduces interactions between selection decisions that preclude exact ILP solution. We solve this formulation by exhaustive enumeration, which remains tractable when $k$ and the number of candidates are small.