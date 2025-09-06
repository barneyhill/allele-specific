#!/usr/bin/env python3
"""
Optimal ASO selection using a greedy algorithm that correctly models
haplotype-specific expected population coverage, restricted to heterozygous sites,
with added summaries for superpopulations and heterozygous proportion tables.
Now includes coverage plots for different k values.
"""

import cyvcf2
import pyensembl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from ortools.sat.python import cp_model

def load_superpopulation_data(vcf_samples):
    """
    Loads superpopulation data and maps it to the samples in the VCF.
    """
    # Load the metadata file
    info_file = 'hgdp1kgp_phased_haplotypes/hgdp_1kg_sample_info.unrelateds.pca_outliers_removed.with_project.tsv'
    df = pd.read_csv(info_file, sep='\t')
    sample_to_pop = dict(zip(df.Sample, df.SuperPop))

    # Map VCF samples to superpopulations
    superpops = [sample_to_pop.get(s) for s in vcf_samples]
    
    # Get a list of unique superpopulations present in the VCF
    unique_pops = list(set(p for p in superpops if p is not None))
    
    # Create boolean masks for each superpopulation
    pop_masks = {pop: np.array([p == pop for p in superpops]) for pop in unique_pops}
    
    return superpops, unique_pops, pop_masks

def get_variants_and_haplotype_matrices(gene_name='SCN2A'):
    """
    Extracts variants and builds two matrices representing the alleles
    on each of the two phased haplotypes for every individual.
    """
    # Get gene coordinates
    ensembl = pyensembl.EnsemblRelease(release=110, species='human')
    gene = ensembl.genes_by_name(gene_name)[0]
    chrom, start, end = gene.contig, gene.start, gene.end

    # Read BCF for gene region
    bcf_path = f"hgdp1kgp_phased_haplotypes/hgdp1kgp_chr{chrom}.filtered.SNV_INDEL.phased.shapeit5.bcf"
    vcf = cyvcf2.VCF(bcf_path)
    
    vcf_samples = vcf.samples

    variants = []
    hap1_vectors, hap2_vectors = [], []

    for variant in vcf(f"chr{chrom}:{start}-{end}"):
        if not variant.is_indel: continue

        variants.append({
            'pos': variant.POS,
            'ref': variant.REF,
            'alt': variant.ALT[0],
        })
        # 0 = REF, 1 = ALT.
        genotypes = variant.genotype.array()[:, 0:2]
        hap1_vectors.append(genotypes[:, 0] == 1)
        hap2_vectors.append(genotypes[:, 1] == 1)

    # Transpose so that rows are individuals and columns are variants
    hap1_matrix = np.array(hap1_vectors, dtype=bool).T
    hap2_matrix = np.array(hap2_vectors, dtype=bool).T

    return variants, hap1_matrix, hap2_matrix, vcf_samples

def ilp_selection_ortools(variants, hap1_matrix, hap2_matrix, k, pop_masks=None):
    """
    Here we pose the problem of selecting k ASOs to maximize population coverage as an integer linear program (ILP).
    We use Google's OR-Tools CP-SAT solver to solve the ILP.
    """

    n_individuals, n_variants = hap1_matrix.shape
    
    # Only count heterozygous sites for each individual
    het_coverage = []
    for j in range(n_individuals):
        hap1_sites = []
        hap2_sites = []
        for i in range(n_variants):
            if hap1_matrix[j, i] != hap2_matrix[j, i]:
                if hap1_matrix[j, i]:
                    hap1_sites.append(i)
                if hap2_matrix[j, i]:
                    hap2_sites.append(i)
        het_coverage.append((hap1_sites, hap2_sites))
    
    # Create CP-SAT model
    model = cp_model.CpModel()
    
    # Binary variables for ASO selection
    x = [model.NewBoolVar(f'x_{i}') for i in range(n_variants)]
    
    # Binary variables for individual coverage
    y1 = [model.NewBoolVar(f'y1_{j}') for j in range(n_individuals)]
    y2 = [model.NewBoolVar(f'y2_{j}') for j in range(n_individuals)]
    
    # Constraint: select at most k ASOs
    model.Add(sum(x) <= k)
    
    # Coverage constraints
    for j in range(n_individuals):
        hap1_sites, hap2_sites = het_coverage[j]
        
        if hap1_sites:
            # y1[j] = 1 if any of the hap1_sites are selected
            model.AddMaxEquality(y1[j], [x[i] for i in hap1_sites])
        else:
            model.Add(y1[j] == 0)
            
        if hap2_sites:
            model.AddMaxEquality(y2[j], [x[i] for i in hap2_sites])
        else:
            model.Add(y2[j] == 0)
    
    # Objective: maximize coverage (scaled by 2 to avoid fractions)
    model.Maximize(sum(y1) + sum(y2))
    
    # Solve with all cores
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 0  # Use all available cores
    solver.parameters.max_time_in_seconds = 300
    solver.parameters.log_search_progress = False
    
    status = solver.Solve(model)
    
    if status in [cp_model.OPTIMAL, cp_model.FEASIBLE]:
        selected_indices = [i for i in range(n_variants) if solver.Value(x[i])]
        coverage = solver.ObjectiveValue() / 2  # Divide by 2 since we scaled
        print(f"OR-Tools found solution with coverage: {coverage}")
        return [variants[i] for i in selected_indices]
    else:
        print("No solution found")
        return []

def create_heterozygous_table(selected_indices, variants, hap1_matrix, hap2_matrix, k):
    """
    Create a k x k table showing heterozygous proportions.
    Diagonal: proportion heterozygous for each variant individually
    Off-diagonal: proportion heterozygous for both variants in the pair
    """
    n_individuals = hap1_matrix.shape[0]
    
    # Sort variants by genomic position for consistent ordering
    sorted_pairs = [(idx, variants[idx]) for idx in selected_indices]
    sorted_pairs.sort(key=lambda x: x[1]['pos'])
    sorted_indices = [pair[0] for pair in sorted_pairs]
    
    # Create the heterozygous proportion matrix
    het_matrix = np.zeros((k, k))
    
    for i in range(k):
        for j in range(k):
            idx_i = sorted_indices[i]
            idx_j = sorted_indices[j]
            
            if i == j:
                # Diagonal: proportion heterozygous for variant i
                is_het_i = hap1_matrix[:, idx_i] != hap2_matrix[:, idx_i]
                het_matrix[i, j] = np.sum(is_het_i) / n_individuals
            else:
                # Off-diagonal: proportion heterozygous for BOTH variants i and j
                is_het_i = hap1_matrix[:, idx_i] != hap2_matrix[:, idx_i]
                is_het_j = hap1_matrix[:, idx_j] != hap2_matrix[:, idx_j]
                both_het = is_het_i & is_het_j
                het_matrix[i, j] = np.sum(both_het) / n_individuals
    
    # Create row/column labels
    labels = []
    for idx in sorted_indices:
        v = variants[idx]
        labels.append(f"chr2:{v['pos']:,}")
    
    # Create DataFrame
    het_df = pd.DataFrame(het_matrix, index=labels, columns=labels)
    
    return het_df

def create_coverage_plots(coverage_data, unique_pops):
    """
    Create and save two scatter plots:
    1. Overall coverage across all populations
    2. Individual lines for each subpopulation
    """
    k_values = list(coverage_data.keys())
    
    # Set up matplotlib style
    plt.style.use('default')
    fig_size = (10, 6)
    
    # Plot 1: Overall coverage
    plt.figure(figsize=fig_size)
    overall_proportions = [coverage_data[k]['overall_proportion'] for k in k_values]
    plt.scatter(k_values, overall_proportions, s=50, alpha=0.7, color='navy')
    plt.plot(k_values, overall_proportions, '-', linewidth=2, color='navy', label='Overall Coverage')
    
    plt.xlabel('Number of ASOs (k)', fontsize=12)
    plt.ylabel('Proportion of Individuals Covered', fontsize=12)
    plt.title('ASO Coverage vs Number of Selected Variants\n(Overall Population)', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Format y-axis as percentages
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: '{:.1%}'.format(y)))
    
    # Set reasonable axis limits
    plt.xlim(0.5, max(k_values) + 0.5)
    plt.ylim(0, max(overall_proportions) * 1.1)
    
    plt.tight_layout()
    plt.savefig('aso_coverage_overall.png', dpi=300, bbox_inches='tight')
    plt.savefig('aso_coverage_overall.pdf', bbox_inches='tight')
    plt.close()
    
    # Plot 2: Individual subpopulation lines
    plt.figure(figsize=fig_size)
    
    # Define colors for different populations
    colors = plt.cm.Set1(np.linspace(0, 1, len(unique_pops)))
    
    for i, pop in enumerate(unique_pops):
        pop_proportions = [coverage_data[k]['pop_proportions'][pop] for k in k_values]
        plt.scatter(k_values, pop_proportions, s=40, alpha=0.7, color=colors[i])
        plt.plot(k_values, pop_proportions, '-', linewidth=2, color=colors[i], label=pop)
    
    plt.xlabel('Number of ASOs (k)', fontsize=12)
    plt.ylabel('Proportion of Individuals Covered', fontsize=12)
    plt.title('ASO Coverage vs Number of Selected Variants\n(By Superpopulation)', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Format y-axis as percentages
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: '{:.1%}'.format(y)))
    
    # Set reasonable axis limits
    plt.xlim(0.5, max(k_values) + 0.5)
    all_proportions = []
    for k in k_values:
        all_proportions.extend(coverage_data[k]['pop_proportions'].values())
    plt.ylim(0, max(all_proportions) * 1.1)
    
    plt.tight_layout()
    plt.savefig('aso_coverage_by_population.png', dpi=300, bbox_inches='tight')
    plt.savefig('aso_coverage_by_population.pdf', bbox_inches='tight')
    plt.close()
    
    print("Coverage plots saved:")
    print("- aso_coverage_overall.png/pdf (Overall population)")
    print("- aso_coverage_by_population.png/pdf (By superpopulation)")

def main():
    """Main function to run ILP ASO selection for k=1 to 5 and display variant coverage table at each k."""
    print("Extracting SCN2A variants from phased BCF...")
    variants, hap1_matrix, hap2_matrix, vcf_samples = get_variants_and_haplotype_matrices('SCN2A')
    n_total = hap1_matrix.shape[0]
    print(f"Found {len(variants)} indel sites across {n_total} individuals\n")

    print("Loading and mapping superpopulation data...")
    _, unique_pops, pop_masks = load_superpopulation_data(vcf_samples)
    pop_counts = {pop: np.sum(mask) for pop, mask in pop_masks.items()}
    print(f"Found {len(unique_pops)} superpopulations: {', '.join(unique_pops)}\n")

    # Helper function to calculate individual variant coverage
    def calculate_variant_coverage(variant_idx, hap1_matrix, hap2_matrix, pop_masks, pop_counts):
        """Calculate coverage provided by a single variant across populations."""
        n_individuals = hap1_matrix.shape[0]
        
        # Identify heterozygous individuals for this variant (only these contribute to ASO coverage)
        is_het = hap1_matrix[:, variant_idx] != hap2_matrix[:, variant_idx]
        hap1_covered = hap1_matrix[:, variant_idx] & is_het
        hap2_covered = hap2_matrix[:, variant_idx] & is_het
        
        total_coverage = 0.5 * (np.sum(hap1_covered) + np.sum(hap2_covered))
        total_proportion = total_coverage / n_individuals
        
        # Calculate per-population coverage as proportions
        pop_proportions = {}
        for pop, mask in pop_masks.items():
            pop_hap1 = np.sum(hap1_covered & mask)
            pop_hap2 = np.sum(hap2_covered & mask)
            pop_coverage = 0.5 * (pop_hap1 + pop_hap2)
            pop_proportions[pop] = pop_coverage / pop_counts[pop] if pop_counts[pop] > 0 else 0.0
        
        return total_proportion, pop_proportions

    # Store coverage data for plotting
    coverage_data = {}

    # Run ILP for different k values and display table after each solve
    for k in range(1, 10):
        print(f"Solving ILP for k={k}...")
        selected = ilp_selection_ortools(variants, hap1_matrix, hap2_matrix, k, pop_masks)
        selected_indices = [variants.index(v) for v in selected]
        
        # Calculate total coverage for this k
        n_individuals = hap1_matrix.shape[0]
        hap1_covered = np.zeros(n_individuals, dtype=bool)
        hap2_covered = np.zeros(n_individuals, dtype=bool)
        
        for idx in selected_indices:
            is_het = hap1_matrix[:, idx] != hap2_matrix[:, idx]
            hap1_covered |= (hap1_matrix[:, idx] & is_het)
            hap2_covered |= (hap2_matrix[:, idx] & is_het)
        
        total_coverage = 0.5 * (np.sum(hap1_covered) + np.sum(hap2_covered))
        overall_proportion = total_coverage / n_total
        
        # Calculate per-population coverage proportions for plotting
        pop_proportions = {}
        for pop in unique_pops:
            mask = pop_masks[pop]
            pop_hap1 = np.sum(hap1_covered & mask)
            pop_hap2 = np.sum(hap2_covered & mask)
            pop_coverage = 0.5 * (pop_hap1 + pop_hap2)
            pop_total = pop_counts[pop]
            pop_proportions[pop] = pop_coverage / pop_total if pop_total > 0 else 0
        
        # Store coverage data for plotting
        coverage_data[k] = {
            'overall_proportion': overall_proportion,
            'pop_proportions': pop_proportions
        }
        
        print(f"OR-Tools found solution with coverage: {total_coverage}")
        
        # Create and display heterozygous proportion table
        if len(selected_indices) > 0:
            print(f"\n=== HETEROZYGOUS PROPORTION TABLE AT k={k} ===")
            het_table = create_heterozygous_table(selected_indices, variants, hap1_matrix, hap2_matrix, k)
            
            # Format as percentages
            het_table_formatted = het_table.applymap(lambda x: f"{x:.1%}")
            print(het_table_formatted.to_string())
            print()
        
        # Create DataFrame for this k's selected variants
        variant_data = []
        
        # Sort variants by genomic position
        sorted_indices = sorted(selected_indices, key=lambda idx: variants[idx]['pos'])
        
        for variant_idx in sorted_indices:
            v = variants[variant_idx]
            total_prop, pop_props = calculate_variant_coverage(variant_idx, hap1_matrix, hap2_matrix, pop_masks, pop_counts)
            
            row_data = {
                'Variant': f"chr2:{v['pos']:,}",
                'Position': v['pos'],
                'Ref/Alt': f"{v['ref']}/{v['alt']}"
            }
            
            # Add population proportions
            for pop in unique_pops:
                row_data[pop] = pop_props[pop]
            
            row_data['Total'] = total_prop
            variant_data.append(row_data)
        
        # Create and display DataFrame
        if variant_data:
            df = pd.DataFrame(variant_data)
            
            print(f"=== SELECTED VARIANTS AT k={k} ===")
            
            # Format the DataFrame for display
            pd.set_option('display.max_columns', None)
            pd.set_option('display.width', None)
            pd.set_option('display.max_colwidth', None)
            
            # Format population columns and Total to show as percentages
            pop_cols = unique_pops + ['Total']
            for col in pop_cols:
                if col in df.columns:
                    df[col] = df[col].apply(lambda x: f"{x:.1%}")
            
            print(df.to_string(index=False, float_format='%.1f'))
            print()
        
        # Calculate and display cumulative coverage summary
        print(f"=== CUMULATIVE COVERAGE SUMMARY AT k={k} ===")
        summary_data = {'Population': [], 'Coverage': [], 'Total': [], 'Proportion': []}
        
        for pop in unique_pops:
            mask = pop_masks[pop]
            pop_hap1 = np.sum(hap1_covered & mask)
            pop_hap2 = np.sum(hap2_covered & mask)
            pop_coverage = 0.5 * (pop_hap1 + pop_hap2)
            pop_total = pop_counts[pop]
            pop_prop = pop_coverage / pop_total if pop_total > 0 else 0
            
            summary_data['Population'].append(pop)
            summary_data['Coverage'].append(f"{pop_coverage:.1f}")
            summary_data['Total'].append(str(pop_total))
            summary_data['Proportion'].append(f"{pop_prop:.1%}")
        
        # Add overall summary
        summary_data['Population'].append('OVERALL')
        summary_data['Coverage'].append(f"{total_coverage:.1f}")
        summary_data['Total'].append(str(n_total))
        summary_data['Proportion'].append(f"{overall_proportion:.1%}")
        
        summary_df = pd.DataFrame(summary_data)
        print(summary_df.to_string(index=False))
        print("=" * 60)
        print()

    # Create and save the coverage plots
    print("\nCreating coverage plots...")
    create_coverage_plots(coverage_data, unique_pops)

if __name__ == "__main__":
    main()