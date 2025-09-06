## Optimization Framework for Maximizing Expected Population Coverage

### **Objective**

For a single target gene, identify the optimal set of `k` ASOs (Antisense Oligonucleotides) that maximizes the **expected number of individuals covered** in a phased haplotype population, where coverage is restricted to heterozygous variants to ensure allele-specific targeting.

### **The Haplotype-Specific Model**

An individual is considered covered if at least one of the selected ASOs targets a variant where:
1. The individual is **heterozygous** for that variant (different alleles on the two haplotypes)
2. The variant appears on at least one of their haplotypes

Key principles:
- Only heterozygous sites contribute to ASO coverage (homozygous sites cannot provide allele-specific targeting)
- Each heterozygous individual has a 50% chance of coverage per haplotype
- Total expected coverage = 0.5 × (individuals with variant on haplotype 1) + 0.5 × (individuals with variant on haplotype 2)

### **The Solution: Integer Linear Programming (ILP)**

The implementation uses Google's OR-Tools CP-SAT solver to find the optimal ASO set through integer linear programming, which guarantees finding the optimal solution (unlike greedy approaches).

#### 1. Input Data

- **Phased Haplotype Data**: BCF files containing phased genotypes from HGDP+1KGP dataset
- **Gene Coordinates**: Retrieved using PyEnsembl (default: SCN2A gene)
- **Population Metadata**: Superpopulation assignments for stratified analysis
- **Variant Information**: Indel positions, reference and alternative alleles

#### 2. Data Structures

- **Haplotype Matrices (`hap1_matrix`, `hap2_matrix`)**: Binary matrices where:
  - Rows represent individuals
  - Columns represent variant positions
  - Values indicate presence (1) or absence (0) of the alternative allele
- **Heterozygous Coverage Lists**: For each individual, tracks which variants are heterozygous and on which haplotype

#### 3. ILP Formulation

**Decision Variables:**
- `x_i`: Binary variable indicating whether ASO `i` is selected
- `y1_j`, `y2_j`: Binary variables indicating whether individual `j` is covered through haplotype 1 or 2

**Constraints:**
- `Σ x_i ≤ k`: Select at most k ASOs
- Coverage constraints: An individual is covered on a haplotype if at least one selected ASO targets a heterozygous variant on that haplotype

**Objective:**
- Maximize: `Σ (y1_j + y2_j)` (scaled expected coverage across both haplotypes)

#### 4. Implementation Features

**Core Analysis:**
- Iterative solving for k=1 to 9 ASOs
- Optimal selection using OR-Tools CP-SAT solver with multi-core support
- Coverage calculation restricted to heterozygous sites only

**Visualization & Reporting:**
- **Heterozygous Proportion Tables**: k×k matrices showing:
  - Diagonal: proportion of individuals heterozygous for each variant
  - Off-diagonal: proportion heterozygous for both variants in a pair
- **Coverage Plots**: 
  - Overall population coverage vs. number of ASOs
  - Stratified coverage by superpopulation
- **Variant Tables**: Selected variants with per-population coverage statistics
- **Cumulative Coverage Summaries**: Total and proportional coverage by population

#### 5. Output

For each value of k, the program produces:
1. Optimal set of k ASO target positions
2. Heterozygous proportion matrix for selected variants
3. Per-variant coverage breakdown by superpopulation
4. Cumulative coverage statistics (total individuals covered and proportions)
5. Publication-quality plots showing coverage trends

### **Usage**

```bash
python main.py
```

The script will:
1. Extract variants from the phased BCF file for the target gene
2. Load superpopulation metadata
3. Solve the ILP optimization for k=1 to 9
4. Display detailed tables and statistics for each k
5. Generate coverage plots saved as PNG and PDF files

### **Dependencies**

- `cyvcf2`: VCF/BCF file parsing
- `pyensembl`: Gene coordinate retrieval
- `numpy`: Matrix operations
- `pandas`: Data manipulation and display
- `matplotlib`: Visualization
- `ortools`: ILP solver (CP-SAT)