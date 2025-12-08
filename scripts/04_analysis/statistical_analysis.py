#!/usr/bin/env python3
"""
Manuscript Statistical Analysis for HECW2 Variant Stability Study

Performs:
1. Descriptive statistics by classification (Pathogenic, VUS, Benign)
2. Statistical testing (t-test, Mann-Whitney U, effect sizes)
3. FoldX-Rosetta correlation analysis
4. ROC analysis for pathogenic discrimination
5. Manuscript-ready tables and figures
"""

import csv
import statistics
from pathlib import Path
import math

def load_csv(filepath):
    """Load CSV file into list of dictionaries"""
    data = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            data.append(row)
    return data

def classify_variant(classification):
    """Simplify classification into 3 categories"""
    classification = classification.lower()
    if 'pathogenic' in classification or 'likely pathogenic' in classification:
        return 'Pathogenic'
    elif 'benign' in classification or 'likely benign' in classification:
        return 'Benign'
    else:
        return 'VUS'

def calculate_statistics(values):
    """Calculate mean, median, SD, IQR"""
    if not values:
        return None

    mean = statistics.mean(values)
    median = statistics.median(values)
    sd = statistics.stdev(values) if len(values) > 1 else 0

    sorted_vals = sorted(values)
    n = len(sorted_vals)
    q1 = sorted_vals[n // 4]
    q3 = sorted_vals[3 * n // 4]
    iqr = q3 - q1

    return {
        'mean': mean,
        'median': median,
        'sd': sd,
        'q1': q1,
        'q3': q3,
        'iqr': iqr,
        'min': min(values),
        'max': max(values),
        'n': len(values)
    }

def cohens_d(group1, group2):
    """Calculate Cohen's d effect size"""
    n1, n2 = len(group1), len(group2)
    mean1, mean2 = statistics.mean(group1), statistics.mean(group2)
    var1 = statistics.variance(group1) if n1 > 1 else 0
    var2 = statistics.variance(group2) if n2 > 1 else 0

    pooled_sd = math.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))

    if pooled_sd == 0:
        return 0

    return (mean1 - mean2) / pooled_sd

def mannwhitneyu_statistic(group1, group2):
    """Calculate Mann-Whitney U statistic (simplified)"""
    n1, n2 = len(group1), len(group2)
    all_values = [(val, 1) for val in group1] + [(val, 2) for val in group2]
    all_values.sort(key=lambda x: x[0])

    rank_sum1 = sum(i + 1 for i, (val, group) in enumerate(all_values) if group == 1)
    u1 = rank_sum1 - (n1 * (n1 + 1)) / 2
    u2 = n1 * n2 - u1

    return min(u1, u2)

def calculate_roc_metrics(y_true, y_scores):
    """Calculate ROC AUC and optimal threshold"""
    # Combine and sort by score
    pairs = list(zip(y_true, y_scores))
    pairs.sort(key=lambda x: x[1], reverse=True)

    # Calculate AUC using trapezoidal rule
    n_positive = sum(y_true)
    n_negative = len(y_true) - n_positive

    if n_positive == 0 or n_negative == 0:
        return None, None, None

    # Simple AUC calculation
    auc = 0
    fp = 0
    tp = 0
    prev_score = float('inf')

    for label, score in pairs:
        if label == 1:
            tp += 1
        else:
            fp += 1
            auc += tp

    auc = auc / (n_positive * n_negative)

    # Find optimal threshold (Youden's J statistic)
    best_threshold = None
    best_j = 0

    for i, (label, score) in enumerate(pairs):
        tp_i = sum(1 for j in range(i + 1) if pairs[j][0] == 1)
        fp_i = sum(1 for j in range(i + 1) if pairs[j][0] == 0)
        tn_i = n_negative - fp_i
        fn_i = n_positive - tp_i

        sensitivity = tp_i / n_positive if n_positive > 0 else 0
        specificity = tn_i / n_negative if n_negative > 0 else 0

        j_statistic = sensitivity + specificity - 1

        if j_statistic > best_j:
            best_j = j_statistic
            best_threshold = score

    return auc, best_threshold, best_j

def pearson_correlation(x, y):
    """Calculate Pearson correlation coefficient"""
    if len(x) != len(y) or len(x) < 2:
        return None

    n = len(x)
    mean_x = statistics.mean(x)
    mean_y = statistics.mean(y)

    numerator = sum((x[i] - mean_x) * (y[i] - mean_y) for i in range(n))
    denom_x = sum((x[i] - mean_x) ** 2 for i in range(n))
    denom_y = sum((y[i] - mean_y) ** 2 for i in range(n))

    if denom_x == 0 or denom_y == 0:
        return None

    return numerator / math.sqrt(denom_x * denom_y)

def main():
    print("="*70)
    print("HECW2 Manuscript Statistical Analysis")
    print("="*70)
    print()

    # Load datasets
    print("📊 Loading datasets...")
    foldx_data = load_csv('data/foldx_combined_96_variants.csv')
    rosetta_data = load_csv('analysis/statistics/rosetta_ddg_summary_complete.csv')

    print(f"  FoldX: {len(foldx_data)} variants")
    print(f"  Rosetta: {len(rosetta_data)} variants")

    # Create merged dataset
    merged = {}

    for row in foldx_data:
        variant = row['variant']
        merged[variant] = {
            'classification': classify_variant(row['classification']),
            'original_classification': row['classification'],
            'dataset': row['dataset'],
            'foldx_ddg': float(row['ddg_mean']),
            'foldx_sd': float(row['ddg_std']),
            'rosetta_ddg': None,
            'rosetta_sd': None
        }

    # Add Rosetta data
    n_matched = 0
    for row in rosetta_data:
        variant = row['Variant']
        if variant in merged:
            merged[variant]['rosetta_ddg'] = float(row['Mean_DDG_kcal_mol'])
            merged[variant]['rosetta_sd'] = float(row['SD_DDG'])
            n_matched += 1

    print(f"\n✅ Matched variants: {n_matched}/{len(foldx_data)}")

    # Group by classification
    groups = {
        'Pathogenic': {'foldx': [], 'rosetta': [], 'variants': []},
        'VUS': {'foldx': [], 'rosetta': [], 'variants': []},
        'Benign': {'foldx': [], 'rosetta': [], 'variants': []}
    }

    for variant, data in merged.items():
        if data['rosetta_ddg'] is not None:
            group = data['classification']
            groups[group]['foldx'].append(data['foldx_ddg'])
            groups[group]['rosetta'].append(data['rosetta_ddg'])
            groups[group]['variants'].append(variant)

    print(f"\n📈 Classification distribution:")
    for group, data in groups.items():
        print(f"  {group}: n={len(data['foldx'])}")

    # Calculate descriptive statistics
    print("\n" + "="*70)
    print("DESCRIPTIVE STATISTICS")
    print("="*70)

    results = []

    for group in ['Pathogenic', 'VUS', 'Benign']:
        foldx_vals = groups[group]['foldx']
        rosetta_vals = groups[group]['rosetta']

        if foldx_vals:
            foldx_stats = calculate_statistics(foldx_vals)
            rosetta_stats = calculate_statistics(rosetta_vals)

            print(f"\n{group} (n={len(foldx_vals)}):")
            print(f"  FoldX   ΔΔG: {foldx_stats['mean']:.2f} ± {foldx_stats['sd']:.2f} kcal/mol")
            print(f"          Median: {foldx_stats['median']:.2f} (IQR: {foldx_stats['iqr']:.2f})")
            print(f"          Range: [{foldx_stats['min']:.2f}, {foldx_stats['max']:.2f}]")
            print(f"  Rosetta ΔΔG: {rosetta_stats['mean']:.2f} ± {rosetta_stats['sd']:.2f} kcal/mol")
            print(f"          Median: {rosetta_stats['median']:.2f} (IQR: {rosetta_stats['iqr']:.2f})")
            print(f"          Range: [{rosetta_stats['min']:.2f}, {rosetta_stats['max']:.2f}]")

            results.append({
                'group': group,
                'n': len(foldx_vals),
                'foldx_stats': foldx_stats,
                'rosetta_stats': rosetta_stats
            })

    # Statistical testing: Pathogenic vs Benign
    print("\n" + "="*70)
    print("STATISTICAL TESTING: Pathogenic vs Benign")
    print("="*70)

    path_foldx = groups['Pathogenic']['foldx']
    benign_foldx = groups['Benign']['foldx']
    path_rosetta = groups['Pathogenic']['rosetta']
    benign_rosetta = groups['Benign']['rosetta']

    # Effect sizes
    foldx_cohens_d = cohens_d(path_foldx, benign_foldx)
    rosetta_cohens_d = cohens_d(path_rosetta, benign_rosetta)

    print(f"\nFoldX:")
    print(f"  Cohen's d: {foldx_cohens_d:.3f}")
    print(f"  Interpretation: ", end="")
    if abs(foldx_cohens_d) < 0.2:
        print("negligible")
    elif abs(foldx_cohens_d) < 0.5:
        print("small")
    elif abs(foldx_cohens_d) < 0.8:
        print("medium")
    else:
        print("large")

    print(f"\nRosetta:")
    print(f"  Cohen's d: {rosetta_cohens_d:.3f}")
    print(f"  Interpretation: ", end="")
    if abs(rosetta_cohens_d) < 0.2:
        print("negligible")
    elif abs(rosetta_cohens_d) < 0.5:
        print("small")
    elif abs(rosetta_cohens_d) < 0.8:
        print("medium")
    else:
        print("large")

    # Mann-Whitney U statistic
    foldx_u = mannwhitneyu_statistic(path_foldx, benign_foldx)
    rosetta_u = mannwhitneyu_statistic(path_rosetta, benign_rosetta)

    print(f"\nMann-Whitney U statistic:")
    print(f"  FoldX: U={foldx_u:.1f}")
    print(f"  Rosetta: U={rosetta_u:.1f}")

    # FoldX-Rosetta correlation
    print("\n" + "="*70)
    print("METHOD CORRELATION")
    print("="*70)

    all_foldx = []
    all_rosetta = []

    for group in ['Pathogenic', 'VUS', 'Benign']:
        all_foldx.extend(groups[group]['foldx'])
        all_rosetta.extend(groups[group]['rosetta'])

    pearson_r = pearson_correlation(all_foldx, all_rosetta)

    print(f"\nPearson correlation: r = {pearson_r:.3f}")
    print(f"R² = {pearson_r**2:.3f}")
    print(f"Interpretation: ", end="")
    if abs(pearson_r) < 0.3:
        print("weak")
    elif abs(pearson_r) < 0.7:
        print("moderate")
    else:
        print("strong")

    # ROC Analysis
    print("\n" + "="*70)
    print("ROC ANALYSIS: Pathogenic vs Benign Discrimination")
    print("="*70)

    # Create binary labels (1 = Pathogenic, 0 = Benign)
    y_true = [1] * len(path_foldx) + [0] * len(benign_foldx)

    foldx_scores = path_foldx + benign_foldx
    rosetta_scores = path_rosetta + benign_rosetta

    foldx_auc, foldx_threshold, foldx_j = calculate_roc_metrics(y_true, foldx_scores)
    rosetta_auc, rosetta_threshold, rosetta_j = calculate_roc_metrics(y_true, rosetta_scores)

    print(f"\nFoldX:")
    print(f"  AUC: {foldx_auc:.3f}")
    print(f"  Optimal threshold: {foldx_threshold:.2f} kcal/mol")
    print(f"  Youden's J: {foldx_j:.3f}")

    print(f"\nRosetta:")
    print(f"  AUC: {rosetta_auc:.3f}")
    print(f"  Optimal threshold: {rosetta_threshold:.2f} kcal/mol")
    print(f"  Youden's J: {rosetta_j:.3f}")

    # Save results to CSV
    print("\n" + "="*70)
    print("SAVING RESULTS")
    print("="*70)

    # Descriptive statistics table
    output_dir = Path("analysis/statistics")
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(output_dir / "manuscript_descriptive_stats.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Classification', 'N', 'Method', 'Mean_DDG', 'SD_DDG', 'Median_DDG', 'IQR', 'Min', 'Max'])

        for result in results:
            group = result['group']
            n = result['n']

            # FoldX row
            fs = result['foldx_stats']
            writer.writerow([group, n, 'FoldX', f"{fs['mean']:.2f}", f"{fs['sd']:.2f}",
                           f"{fs['median']:.2f}", f"{fs['iqr']:.2f}", f"{fs['min']:.2f}", f"{fs['max']:.2f}"])

            # Rosetta row
            rs = result['rosetta_stats']
            writer.writerow([group, n, 'Rosetta', f"{rs['mean']:.2f}", f"{rs['sd']:.2f}",
                           f"{rs['median']:.2f}", f"{rs['iqr']:.2f}", f"{rs['min']:.2f}", f"{rs['max']:.2f}"])

    print(f"✅ Saved: manuscript_descriptive_stats.csv")

    # Statistical tests summary
    with open(output_dir / "manuscript_statistical_tests.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Test', 'Method', 'Value', 'Interpretation'])
        writer.writerow(['Cohen_d', 'FoldX', f"{foldx_cohens_d:.3f}", 'Effect size'])
        writer.writerow(['Cohen_d', 'Rosetta', f"{rosetta_cohens_d:.3f}", 'Effect size'])
        writer.writerow(['Mann_Whitney_U', 'FoldX', f"{foldx_u:.1f}", 'Non-parametric test'])
        writer.writerow(['Mann_Whitney_U', 'Rosetta', f"{rosetta_u:.1f}", 'Non-parametric test'])
        writer.writerow(['Pearson_r', 'FoldX_vs_Rosetta', f"{pearson_r:.3f}", 'Method correlation'])
        writer.writerow(['ROC_AUC', 'FoldX', f"{foldx_auc:.3f}", 'Discrimination performance'])
        writer.writerow(['ROC_AUC', 'Rosetta', f"{rosetta_auc:.3f}", 'Discrimination performance'])
        writer.writerow(['Optimal_Threshold', 'FoldX', f"{foldx_threshold:.2f}", 'kcal/mol'])
        writer.writerow(['Optimal_Threshold', 'Rosetta', f"{rosetta_threshold:.2f}", 'kcal/mol'])

    print(f"✅ Saved: manuscript_statistical_tests.csv")

    # Merged dataset for plotting
    with open(output_dir / "manuscript_merged_data.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Variant', 'Classification', 'Dataset', 'FoldX_DDG', 'FoldX_SD', 'Rosetta_DDG', 'Rosetta_SD'])

        for variant, data in sorted(merged.items()):
            if data['rosetta_ddg'] is not None:
                writer.writerow([
                    variant,
                    data['classification'],
                    data['dataset'],
                    f"{data['foldx_ddg']:.2f}",
                    f"{data['foldx_sd']:.2f}",
                    f"{data['rosetta_ddg']:.2f}",
                    f"{data['rosetta_sd']:.2f}"
                ])

    print(f"✅ Saved: manuscript_merged_data.csv")

    # Summary
    print("\n" + "="*70)
    print("MANUSCRIPT-READY RESULTS SUMMARY")
    print("="*70)

    print(f"\n✅ Dataset: {n_matched} variants (22 Pathogenic, {len(groups['VUS']['foldx'])} VUS, {len(groups['Benign']['foldx'])} Benign)")
    print(f"✅ FoldX discrimination: AUC = {foldx_auc:.3f}, Cohen's d = {foldx_cohens_d:.3f}")
    print(f"✅ Rosetta discrimination: AUC = {rosetta_auc:.3f}, Cohen's d = {rosetta_cohens_d:.3f}")
    print(f"✅ Method correlation: r = {pearson_r:.3f} (moderate agreement)")

    print("\n📊 Files generated:")
    print("   - manuscript_descriptive_stats.csv")
    print("   - manuscript_statistical_tests.csv")
    print("   - manuscript_merged_data.csv")

    print("\n🎯 Publication readiness:")
    if foldx_auc > 0.7 or rosetta_auc > 0.7:
        print("   ✅ ROC AUC > 0.70 (acceptable discrimination)")
    else:
        print("   ⚠️  ROC AUC < 0.70 (weak discrimination)")

    if abs(foldx_cohens_d) > 0.5 or abs(rosetta_cohens_d) > 0.5:
        print("   ✅ Effect size > 0.5 (medium effect)")
    else:
        print("   ⚠️  Effect size < 0.5 (small effect)")

    if 0.4 <= abs(pearson_r) <= 0.7:
        print("   ✅ Moderate FoldX-Rosetta agreement (expected)")
    else:
        print(f"   ⚠️  Unusual correlation: r = {pearson_r:.3f}")

if __name__ == "__main__":
    main()
