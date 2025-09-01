## Optimization Framework for Maximizing Expected Population Coverage

### **Objective**

For a single target gene, identify the optimal minimal set of `K` ASOs that maximizes the **expected number of individuals covered** in a population, assuming a 50% probability that a targeted heterozygous variant is on the disease-associated haplotype.

### **The Probabilistic Model**

An individual is considered successfully covered if at least one of the selected ASOs targets a variant that is both heterozygous in that individual and resides on the disease haplotype.

*   The probability that a single suitable ASO fails (i.e., is on the healthy haplotype) is `0.5`.
*   The probability that an individual `i` is **not covered** by a selected set of ASOs `S` is the probability that *all* suitable ASOs in that set fail simultaneously.
*   The probability that an individual `i` **is covered** is `1 - P(not covered)`.

If an individual `i` is heterozygous for `c` ASOs within the selected set `S`, their probability of being covered is:
`P(Coveredᵢ) = 1 - (0.5)ᶜ`

The total expected number of individuals covered for a set `S` is the sum of these probabilities across the entire population:
`ExpectedCoverage(S) = Σ (1 - 0.5ᶜᵢ)`

### **The Solution: A Greedy Algorithm**

A greedy algorithm provides an elegant and computationally efficient solution by iteratively building the optimal set. At each step, it selects the single best ASO that adds the most to the total expected coverage.

#### 1. Given Inputs

*   **Set of Individuals (`I`)**: A set representing all `n` individuals in the cohort.
*   **Set of Candidate ASOs (`A`)**: A set of all `m` potential ASOs for the gene.
*   **Heterozygosity Matrix (`H_ij`)**: A binary matrix derived from WGS data, indicating which individuals are heterozygous for the variants targeted by each ASO.
    *   `H_ij = 1` if individual `i` is heterozygous for the variant targeted by ASO `j`.
    *   `H_ij = 0` otherwise.

#### 2. The Algorithm

1.  **Initialization**:
    *   Start with an empty set of selected ASOs: `S = {}`.
    *   Create a list of remaining candidate ASOs: `Candidates = A`.
    *   Initialize the current total expected coverage: `TotalExpectedCoverage = 0`.

2.  **Iterative Selection**:
    *   Repeat `K` times (for a final set of size `K`):
        a.  Set `BestASO = null` and `MaxMarginalGain = -1`.
        b.  **For each `aso_j` in `Candidates`**:
            i.  Calculate the **marginal gain** in expected coverage that `aso_j` would provide if added to the current set `S`.
            ii. The marginal gain is `ExpectedCoverage(S U {aso_j}) - TotalExpectedCoverage`.
            iii. If this gain is greater than `MaxMarginalGain`, update `MaxMarginalGain` and set `BestASO = aso_j`.
        c.  **Add `BestASO` to `S`**.
        d.  **Remove `BestASO` from `Candidates`**.
        e.  **Update `TotalExpectedCoverage`** to `TotalExpectedCoverage + MaxMarginalGain`.
        f.  Store the results for this step (`K`, the new `S`, and the new `TotalExpectedCoverage`).

#### 3. Output

The algorithm produces a performance table showing the optimal set of ASOs and the corresponding expected population coverage for each set size from 1 to `K`.
