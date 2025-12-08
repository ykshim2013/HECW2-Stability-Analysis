#!/usr/bin/env python3
"""
Parse Rosetta cartesian_ddg results and calculate ΔΔG values
Specifically designed for HECW2 project output format
"""

import os
import re
import pandas as pd
import numpy as np
from pathlib import Path

# Directories
BASE_DIR = '/Users/ykshim2025/Desktop/Code2025/HECW2'
ROSETTA_RESULTS_DIR = os.path.join(BASE_DIR, 'analysis', 'rosetta', 'cartesian_ddg', 'results')
OUTPUT_DIR = os.path.join(BASE_DIR, 'analysis', 'statistics')

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Variant classifications
PATHOGENIC = [
    'p.ARG1191GLN', 'p.ARG1330TRP', 'p.PHE1193VAL', 'p.GLU1445GLY',
    'p.GLU1445GLN', 'p.LEU1448SER', 'p.ARG1191TRP'
]
LIKELY_PATHOGENIC = [
    'p.GLU1491GLN', 'p.GLU1564LYS', 'p.GLU1572LYS', 'p.ILE138VAL'
]
VUS = [
    'p.ARG1222LYS', 'p.ASN1244HIS', 'p.ASN1464SER', 'p.ASP1134GLY',
    'p.GLN1232LYS', 'p.GLU1264ASP', 'p.GLU744LYS', 'p.GLY1177ARG',
    'p.GLY1469VAL', 'p.LEU1155PHE', 'p.LEU1496VAL', 'p.MET1234THR',
    'p.PHE1569SER', 'p.PRO1151SER', 'p.PRO1548SER', 'p.SER1551ALA',
    'p.THR69ALA', 'p.VAL1169LEU'
]
BENIGN = ['p.ARG15GLN', 'p.GLU39ALA']

def parse_rosetta_ddg_file(filepath):
    """
    Parse a single Rosetta mutations.ddg file

    Format:
    COMPLEX:   Round1: WT: 15267.858  fa_atr: ... [energy terms]
    COMPLEX:   Round1: MUT_1572ALA: 15272.345  fa_atr: ... [energy terms]
    """

    wt_energies = []
    mut_energies = []

    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith('COMPLEX:'):
                continue

            # Extract round number and type (WT or MUT)
            match_wt = re.search(r'Round(\d+):\s+WT:\s+([\d.]+)', line)
            match_mut = re.search(r'Round(\d+):\s+MUT_\w+:\s+([\d.]+)', line)

            if match_wt:
                round_num = int(match_wt.group(1))
                energy = float(match_wt.group(2))
                wt_energies.append({'round': round_num, 'energy': energy})

            elif match_mut:
                round_num = int(match_mut.group(1))
                energy = float(match_mut.group(2))
                mut_energies.append({'round': round_num, 'energy': energy})

    return wt_energies, mut_energies

def calculate_ddg_from_energies(wt_energies, mut_energies):
    """
    Calculate ΔΔG = <E_mut> - <E_wt>
    Returns mean, SD, and number of measurements
    """

    if not wt_energies or not mut_energies:
        return None, None, 0

    # Convert to arrays
    wt_array = np.array([e['energy'] for e in wt_energies])
    mut_array = np.array([e['energy'] for e in mut_energies])

    # Calculate mean energies
    mean_wt = np.mean(wt_array)
    mean_mut = np.mean(mut_array)

    # ΔΔG = E_mut - E_wt
    ddg = mean_mut - mean_wt

    # Standard deviation (propagation of error)
    # σ_ΔΔG = sqrt(σ_mut^2 + σ_wt^2)
    sd_wt = np.std(wt_array, ddof=1) if len(wt_array) > 1 else 0
    sd_mut = np.std(mut_array, ddof=1) if len(mut_array) > 1 else 0
    sd_ddg = np.sqrt(sd_wt**2 + sd_mut**2)

    n_measurements = min(len(wt_energies), len(mut_energies))

    return ddg, sd_ddg, n_measurements

def get_variant_classification(variant):
    """Classify variant as Pathogenic, Likely_Pathogenic, VUS, or Benign"""
    if variant in PATHOGENIC:
        return 'Pathogenic'
    elif variant in LIKELY_PATHOGENIC:
        return 'Likely_Pathogenic'
    elif variant in VUS:
        return 'VUS'
    elif variant in BENIGN:
        return 'Benign'
    else:
        return 'Unknown'

def main():
    print("=" * 80)
    print("ROSETTA CARTESIAN_DDG RESULT PARSER")
    print("=" * 80)
    print()

    # Find all mutations.ddg files
    results_dirs = [d for d in Path(ROSETTA_RESULTS_DIR).iterdir() if d.is_dir()]

    print(f"Found {len(results_dirs)} variant directories")
    print()

    # Parse all results
    data = []
    failed = []

    for variant_dir in sorted(results_dirs):
        variant_name = variant_dir.name
        ddg_file = variant_dir / 'mutations.ddg'

        if not ddg_file.exists():
            print(f"  ⚠️  Missing: {variant_name}/mutations.ddg")
            failed.append(variant_name)
            continue

        # Extract variant and seed from directory name
        match = re.match(r'(p\.\w+)_seed(\d+)', variant_name)
        if not match:
            print(f"  ⚠️  Cannot parse: {variant_name}")
            failed.append(variant_name)
            continue

        variant = match.group(1)
        seed = int(match.group(2))

        # Parse ddg file
        try:
            wt_energies, mut_energies = parse_rosetta_ddg_file(ddg_file)
            ddg, sd, n = calculate_ddg_from_energies(wt_energies, mut_energies)

            if ddg is None:
                print(f"  ⚠️  No data: {variant_name}")
                failed.append(variant_name)
                continue

            classification = get_variant_classification(variant)

            data.append({
                'Variant': variant,
                'Seed': seed,
                'Method': 'Rosetta',
                'DDG_kcal_mol': ddg,
                'SD_kcal_mol': sd,
                'N_rounds': n,
                'Classification': classification
            })

            print(f"  ✅ {variant_name}: ΔΔG = {ddg:.2f} ± {sd:.2f} kcal/mol ({n} rounds)")

        except Exception as e:
            print(f"  ❌ Error parsing {variant_name}: {e}")
            failed.append(variant_name)

    # Create DataFrame
    df = pd.DataFrame(data)

    print()
    print("=" * 80)
    print("PARSING SUMMARY")
    print("=" * 80)
    print(f"Total directories: {len(results_dirs)}")
    print(f"Successfully parsed: {len(data)}")
    print(f"Failed: {len(failed)}")
    print()

    if failed:
        print("Failed variants:")
        for f in failed:
            print(f"  - {f}")
        print()

    # Save raw data
    csv_file = os.path.join(OUTPUT_DIR, 'rosetta_ddg_all_variants.csv')
    df.to_csv(csv_file, index=False)
    print(f"✅ Raw data saved to: {csv_file}")
    print()

    # Calculate summary statistics by variant (average across seeds)
    summary_data = []

    for variant in df['Variant'].unique():
        variant_data = df[df['Variant'] == variant]

        classification = variant_data['Classification'].iloc[0]

        # Average ΔΔG across seeds
        mean_ddg = variant_data['DDG_kcal_mol'].mean()
        sd_ddg = variant_data['DDG_kcal_mol'].std(ddof=1) if len(variant_data) > 1 else 0
        n_seeds = len(variant_data)

        # Average SD across seeds
        mean_sd = variant_data['SD_kcal_mol'].mean()

        summary_data.append({
            'Variant': variant,
            'Classification': classification,
            'Mean_DDG_kcal_mol': mean_ddg,
            'SD_DDG_kcal_mol': sd_ddg,
            'Mean_Technical_SD': mean_sd,
            'N_Seeds': n_seeds
        })

    summary_df = pd.DataFrame(summary_data)

    # Save summary
    summary_file = os.path.join(OUTPUT_DIR, 'rosetta_ddg_summary.csv')
    summary_df.to_csv(summary_file, index=False)
    print(f"✅ Summary saved to: {summary_file}")
    print()

    # Print classification statistics
    print("=" * 80)
    print("STATISTICS BY CLASSIFICATION")
    print("=" * 80)
    print()

    for classification in ['Pathogenic', 'Likely_Pathogenic', 'VUS', 'Benign']:
        class_data = summary_df[summary_df['Classification'] == classification]

        if len(class_data) == 0:
            continue

        mean = class_data['Mean_DDG_kcal_mol'].mean()
        sd = class_data['Mean_DDG_kcal_mol'].std()
        n = len(class_data)

        print(f"{classification} (n={n}):")
        print(f"  Mean ΔΔG: {mean:.2f} ± {sd:.2f} kcal/mol")
        print(f"  Range: [{class_data['Mean_DDG_kcal_mol'].min():.2f}, {class_data['Mean_DDG_kcal_mol'].max():.2f}]")
        print()

    # Check for destabilizing variants (ΔΔG > 1.0 kcal/mol)
    destabilizing = summary_df[summary_df['Mean_DDG_kcal_mol'] > 1.0].sort_values('Mean_DDG_kcal_mol', ascending=False)

    if len(destabilizing) > 0:
        print("=" * 80)
        print("DESTABILIZING VARIANTS (ΔΔG > 1.0 kcal/mol)")
        print("=" * 80)
        print()

        for idx, row in destabilizing.iterrows():
            print(f"  {row['Variant']:20s} | ΔΔG = {row['Mean_DDG_kcal_mol']:6.2f} ± {row['SD_DDG_kcal_mol']:5.2f} | {row['Classification']}")
        print()

    print("=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print()
    print("Next steps:")
    print("1. Review summary: analysis/statistics/rosetta_ddg_summary.csv")
    print("2. Wait for FoldX results to complete (48-72 hours)")
    print("3. Run full analysis with both methods:")
    print("   python3 scripts/parse_ddg_results.py")
    print()

if __name__ == '__main__':
    main()
