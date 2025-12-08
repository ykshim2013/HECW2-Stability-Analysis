#!/usr/bin/env python3
"""
Parse FoldX BuildModel Results

Extracts ΔΔG values from Dif_*.fxout files for all 3 seeds.
Aggregates results per variant using median ± IQR (robust to outliers).

Usage:
    python3 scripts/parse_foldx_results.py

Output:
    - analysis/statistics/foldx_ddg_all_systems.csv (all measurements)
    - analysis/statistics/foldx_ddg_summary.csv (aggregated per variant)
    - FOLDX_RESULTS_SUMMARY.md (detailed analysis)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import re

# Project paths
PROJECT_DIR = Path("/Users/ykshim2025/Desktop/Code2025/HECW2")
FOLDX_DIR = PROJECT_DIR / "analysis" / "foldx"
STATS_DIR = PROJECT_DIR / "analysis" / "statistics"
STATS_DIR.mkdir(exist_ok=True)

# Variant classifications (from variant_list.csv)
CLASSIFICATIONS = {
    'p.ARG15GLN': 'VUS', 'p.ILE138VAL': 'LP', 'p.GLU39ALA': 'VUS',
    'p.THR69ALA': 'VUS', 'p.GLU744LYS': 'VUS', 'p.ASP1134GLY': 'VUS',
    'p.PRO1151SER': 'P', 'p.LEU1155PHE': 'VUS', 'p.VAL1169LEU': 'VUS',
    'p.GLY1177ARG': 'P', 'p.ARG1191GLN': 'VUS', 'p.PHE1193VAL': 'P',
    'p.GLN1195ARG': 'P', 'p.MET1234THR': 'VUS', 'p.TYR1238ASP': 'P',
    'p.ARG1222LYS': 'VUS', 'p.ILE1332VAL': 'VUS', 'p.ARG1330TRP': 'VUS',
    'p.LEU1410PRO': 'P', 'p.ARG1414GLN': 'VUS', 'p.ARG1432HIS': 'VUS',
    'p.LEU1448SER': 'P', 'p.GLY1469VAL': 'LP', 'p.GLU1491GLN': 'VUS',
    'p.LEU1496VAL': 'VUS', 'p.PRO1548SER': 'VUS', 'p.SER1551ALA': 'VUS',
    'p.GLU1564LYS': 'LP', 'p.PHE1569SER': 'VUS', 'p.GLU1572LYS': 'LP',
    'p.THR1293MET': 'Benign', 'p.SER1447PHE': 'Benign'
}

def parse_foldx_dif_file(filepath):
    """
    Parse FoldX Dif_*.fxout file

    Format:
        PDB file analysed: batch
        Output type: BuildModel
        Pdb	total energy	Backbone Hbond	...
        wt_seed1_1_0.pdb	-0.102171	0.523728	...
        wt_seed1_1_1.pdb	-0.105994	-0.603671	...
        ...

    Returns:
        DataFrame with columns: pdb_file, ddg_kcal_mol, seed, variant_id, run_id
    """
    print(f"Parsing {filepath}")

    # Read file, skip header lines
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Find header line (starts with "Pdb")
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith('Pdb'):
            header_idx = i
            break

    if header_idx is None:
        raise ValueError(f"Could not find header in {filepath}")

    # Extract seed name from filepath
    seed_match = re.search(r'seed(\d)', str(filepath))
    seed = f"seed{seed_match.group(1)}" if seed_match else "unknown"

    # Parse data lines
    data = []
    for line in lines[header_idx + 1:]:
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        parts = line.split('\t')
        if len(parts) < 2:
            continue

        pdb_file = parts[0]
        ddg = float(parts[1])

        # Extract variant_id and run_id from filename
        # Format: wt_seed1_1_0.pdb (variant 1, run 0)
        # Format: WT_wt_seed1_1_0.pdb (wildtype variant 1, run 0)
        match = re.search(r'_(\d+)_(\d+)\.pdb', pdb_file)
        if match:
            variant_id = int(match.group(1))
            run_id = int(match.group(2))
        else:
            print(f"Warning: Could not parse filename: {pdb_file}")
            continue

        data.append({
            'pdb_file': pdb_file,
            'ddg_kcal_mol': ddg,
            'seed': seed,
            'variant_id': variant_id,
            'run_id': run_id
        })

    df = pd.DataFrame(data)
    print(f"  Parsed {len(df)} measurements from {seed}")
    return df

def map_variant_id_to_name():
    """
    Map variant_id (1-31) to variant name based on mutation list order

    Returns:
        dict: {variant_id: variant_name}
    """
    mutation_list = FOLDX_DIR / "individual_list_corrected.txt"

    mapping = {}
    with open(mutation_list, 'r') as f:
        for idx, line in enumerate(f, start=1):
            line = line.strip().rstrip(';')
            if not line:
                continue

            # Format: RA15Q (R at position 15 mutated to Q)
            old_aa_1 = line[0]
            new_aa_1 = line[-1]
            position = line[2:-1]  # Extract position (e.g., "A15" -> "15")
            position = position.lstrip('A')  # Remove chain identifier

            # Convert 1-letter code to 3-letter code
            aa_map = {
                'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
                'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
                'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
                'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
            }

            old_aa = aa_map.get(old_aa_1, old_aa_1)
            new_aa = aa_map.get(new_aa_1, new_aa_1)

            variant_name = f"p.{old_aa}{position}{new_aa}"
            mapping[idx] = variant_name

    return mapping

def aggregate_by_variant(df_all, variant_mapping):
    """
    Aggregate ΔΔG values per variant across all seeds and runs

    Uses median ± IQR for robustness to outliers (FoldX recommendation)

    Returns:
        DataFrame with columns: variant, median_ddg, iqr, classification, etc.
    """
    # Add variant names
    df_all['variant'] = df_all['variant_id'].map(variant_mapping)

    # Calculate statistics per variant
    grouped = df_all.groupby('variant')['ddg_kcal_mol']

    summary = pd.DataFrame({
        'variant': grouped.groups.keys(),
        'median_ddg_kcal_mol': grouped.median(),
        'iqr_ddg': grouped.apply(lambda x: x.quantile(0.75) - x.quantile(0.25)),
        'mean_ddg_kcal_mol': grouped.mean(),
        'sd_ddg': grouped.std(),
        'min_ddg': grouped.min(),
        'max_ddg': grouped.max(),
        'n_measurements': grouped.count(),
        'q1_ddg': grouped.quantile(0.25),
        'q3_ddg': grouped.quantile(0.75)
    }).reset_index(drop=True)

    # Add classifications
    summary['classification'] = summary['variant'].map(CLASSIFICATIONS)

    # Interpret ΔΔG
    def interpret_ddg(ddg):
        if ddg > 1.0:
            return 'Destabilizing'
        elif ddg < -1.0:
            return 'Stabilizing'
        else:
            return 'Neutral'

    summary['interpretation'] = summary['median_ddg_kcal_mol'].apply(interpret_ddg)

    # Sort by median ΔΔG (most destabilizing first)
    summary = summary.sort_values('median_ddg_kcal_mol', ascending=False)

    return summary

def generate_summary_report(df_summary):
    """Generate markdown summary report"""

    report = f"""# FoldX Results Summary

**Date**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}
**Analysis**: FoldX 5.1 BuildModel ΔΔG predictions

---

## Dataset Overview

**Total variants analyzed**: {len(df_summary)}
**Measurements per variant**: 30 (3 seeds × 10 runs)
**Statistical method**: Median ± IQR (robust to outliers)

### Variant Classifications
"""

    for cls in ['P', 'LP', 'VUS', 'Benign']:
        cls_name = {
            'P': 'Pathogenic',
            'LP': 'Likely Pathogenic',
            'VUS': 'Variant of Unknown Significance',
            'Benign': 'Benign/Likely Benign'
        }[cls]
        count = (df_summary['classification'] == cls).sum()
        report += f"- **{cls_name}**: {count} variants\n"

    report += f"\n---\n\n## Top 10 Most Destabilizing Variants\n\n"
    report += "| Rank | Variant | Median ΔΔG | IQR | Classification | Interpretation |\n"
    report += "|------|---------|------------|-----|----------------|----------------|\n"

    for idx, row in df_summary.head(10).iterrows():
        report += f"| {idx+1} | {row['variant']} | {row['median_ddg_kcal_mol']:.2f} | "
        report += f"{row['iqr_ddg']:.2f} | {row['classification']} | {row['interpretation']} |\n"

    report += f"\n---\n\n## Top 10 Most Stabilizing Variants\n\n"
    report += "| Rank | Variant | Median ΔΔG | IQR | Classification | Interpretation |\n"
    report += "|------|---------|------------|-----|----------------|----------------|\n"

    for idx, row in df_summary.tail(10).iterrows():
        report += f"| {idx+1} | {row['variant']} | {row['median_ddg_kcal_mol']:.2f} | "
        report += f"{row['iqr_ddg']:.2f} | {row['classification']} | {row['interpretation']} |\n"

    # Pathogenic vs Benign comparison
    report += f"\n---\n\n## Pathogenic vs Benign Comparison\n\n"

    pathogenic = df_summary[df_summary['classification'].isin(['P', 'LP'])]
    benign = df_summary[df_summary['classification'] == 'Benign']

    report += f"**Pathogenic/Likely Pathogenic** (n={len(pathogenic)}):\n"
    report += f"- Median ΔΔG: {pathogenic['median_ddg_kcal_mol'].median():.2f} kcal/mol\n"
    report += f"- Range: {pathogenic['median_ddg_kcal_mol'].min():.2f} to {pathogenic['median_ddg_kcal_mol'].max():.2f} kcal/mol\n\n"

    report += f"**Benign** (n={len(benign)}):\n"
    report += f"- Median ΔΔG: {benign['median_ddg_kcal_mol'].median():.2f} kcal/mol\n"
    report += f"- Range: {benign['median_ddg_kcal_mol'].min():.2f} to {benign['median_ddg_kcal_mol'].max():.2f} kcal/mol\n\n"

    report += f"\n---\n\n## Quality Control\n\n"
    report += f"**Measurements per variant**: {df_summary['n_measurements'].unique()[0]}\n"
    report += f"**Median IQR**: {df_summary['iqr_ddg'].median():.2f} kcal/mol\n"
    report += f"**Max IQR**: {df_summary['iqr_ddg'].max():.2f} kcal/mol\n"
    report += f"**Variants with IQR > 2.0**: {(df_summary['iqr_ddg'] > 2.0).sum()}\n"

    report += f"\n---\n\n## Next Steps\n\n"
    report += "1. Combine with Rosetta cartesian_ddG results\n"
    report += "2. Calculate FoldX-Rosetta correlation\n"
    report += "3. Generate method comparison plots\n"
    report += "4. Statistical analysis (t-test, ROC curves)\n"
    report += "5. Clinical correlation with disease severity\n"

    return report

def main():
    print("=" * 60)
    print("FoldX Results Parser")
    print("=" * 60)

    # Parse all Dif files
    all_data = []
    for seed in ['seed1', 'seed2', 'seed3']:
        dif_file = FOLDX_DIR / f"{seed}_foldx" / f"Dif_wt_{seed}.fxout"
        if dif_file.exists():
            df_seed = parse_foldx_dif_file(dif_file)
            all_data.append(df_seed)
        else:
            print(f"Warning: {dif_file} not found")

    if not all_data:
        print("Error: No Dif files found!")
        return

    # Combine all seeds
    df_all = pd.concat(all_data, ignore_index=True)
    print(f"\nTotal measurements: {len(df_all)}")

    # Map variant IDs to names
    variant_mapping = map_variant_id_to_name()
    print(f"Mapped {len(variant_mapping)} variants")

    # Aggregate by variant
    print("\nAggregating by variant...")
    df_summary = aggregate_by_variant(df_all, variant_mapping)

    # Save results
    all_csv = STATS_DIR / "foldx_ddg_all_systems.csv"
    summary_csv = STATS_DIR / "foldx_ddg_summary.csv"

    # Add variant names to df_all before saving
    df_all['variant'] = df_all['variant_id'].map(variant_mapping)
    df_all['classification'] = df_all['variant'].map(CLASSIFICATIONS)

    df_all.to_csv(all_csv, index=False)
    df_summary.to_csv(summary_csv, index=False)

    print(f"\n✅ Saved results:")
    print(f"   {all_csv} ({len(df_all)} measurements)")
    print(f"   {summary_csv} ({len(df_summary)} variants)")

    # Generate summary report
    report = generate_summary_report(df_summary)
    report_file = PROJECT_DIR / "FOLDX_RESULTS_SUMMARY.md"
    with open(report_file, 'w') as f:
        f.write(report)
    print(f"   {report_file}")

    print("\n" + "=" * 60)
    print("FoldX Results Parsing Complete!")
    print("=" * 60)

if __name__ == "__main__":
    main()
