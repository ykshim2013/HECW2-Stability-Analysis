#!/usr/bin/env python3
"""
Analyze AlphaFold3 output for HECW2
- Convert CIF to PDB
- Extract pLDDT scores
- Calculate domain-specific confidence
- Generate quality report
"""

import os
import json
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select

# Define HECW2 domains
DOMAINS = {
    'N-terminal': (1, 191),
    'C2': (192, 281),
    'Inter1': (282, 808),
    'WW1': (809, 838),
    'Inter2': (839, 986),
    'WW2': (987, 1016),
    'WW-HECT_linker': (1017, 1267),
    'HECT': (1268, 1571),
    'C-terminal': (1572, 1572)
}

def convert_cif_to_pdb(cif_file, pdb_file):
    """Convert mmCIF to PDB format"""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('HECW2', cif_file)

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)
    print(f"Converted {os.path.basename(cif_file)} → {os.path.basename(pdb_file)}")

def extract_plddt_scores(cif_file):
    """Extract per-residue pLDDT scores from CIF file"""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('HECW2', cif_file)

    plddt_scores = []
    residue_ids = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':  # Standard residue
                    residue_id = residue.id[1]
                    residue_ids.append(residue_id)

                    # pLDDT is stored in B-factor column
                    ca_atoms = [atom for atom in residue if atom.name == 'CA']
                    if ca_atoms:
                        plddt = ca_atoms[0].bfactor
                        plddt_scores.append(plddt)

    return np.array(residue_ids), np.array(plddt_scores)

def calculate_domain_confidence(residue_ids, plddt_scores, domains):
    """Calculate average pLDDT for each domain"""
    domain_stats = {}

    for domain_name, (start, end) in domains.items():
        mask = (residue_ids >= start) & (residue_ids <= end)
        domain_plddt = plddt_scores[mask]

        if len(domain_plddt) > 0:
            domain_stats[domain_name] = {
                'mean': np.mean(domain_plddt),
                'median': np.median(domain_plddt),
                'min': np.min(domain_plddt),
                'max': np.max(domain_plddt),
                'std': np.std(domain_plddt),
                'num_residues': len(domain_plddt)
            }

    return domain_stats

def analyze_model(model_num, base_dir, output_dir):
    """Analyze single AlphaFold model"""
    base_name = 'fold_hecw2_2025_10_10_18_12'
    cif_file = os.path.join(base_dir, f'{base_name}_model_{model_num}.cif')
    pdb_file = os.path.join(output_dir, f'HECW2_model{model_num}.pdb')
    conf_file = os.path.join(base_dir, f'{base_name}_summary_confidences_{model_num}.json')

    # Convert CIF to PDB
    convert_cif_to_pdb(cif_file, pdb_file)

    # Extract pLDDT scores
    residue_ids, plddt_scores = extract_plddt_scores(cif_file)

    # Calculate domain confidence
    domain_stats = calculate_domain_confidence(residue_ids, plddt_scores, DOMAINS)

    # Load summary confidence
    with open(conf_file, 'r') as f:
        summary_conf = json.load(f)

    # Overall statistics
    overall_stats = {
        'mean_plddt': np.mean(plddt_scores),
        'median_plddt': np.median(plddt_scores),
        'num_residues': len(plddt_scores),
        'ptm': summary_conf['ptm'],
        'ranking_score': summary_conf['ranking_score'],
        'fraction_disordered': summary_conf['fraction_disordered'],
        'has_clash': summary_conf['has_clash']
    }

    return {
        'model': model_num,
        'pdb_file': pdb_file,
        'overall': overall_stats,
        'domains': domain_stats,
        'plddt_per_residue': list(zip(residue_ids.tolist(), plddt_scores.tolist()))
    }

def generate_quality_report(results, output_file):
    """Generate comprehensive quality report"""
    with open(output_file, 'w') as f:
        f.write("# HECW2 AlphaFold3 Quality Assessment Report\n")
        f.write(f"# Generated: {os.popen('date').read().strip()}\n")
        f.write("# " + "="*70 + "\n\n")

        # Summary table
        f.write("## Overall Model Quality\n\n")
        f.write("| Model | Mean pLDDT | PTM | Ranking Score | Fraction Disordered | Has Clash |\n")
        f.write("|-------|------------|-----|---------------|--------------------|-----------|\n")

        for result in results:
            f.write(f"| {result['model']} | "
                   f"{result['overall']['mean_plddt']:.2f} | "
                   f"{result['overall']['ptm']:.3f} | "
                   f"{result['overall']['ranking_score']:.3f} | "
                   f"{result['overall']['fraction_disordered']:.2f} | "
                   f"{result['overall']['has_clash']:.1f} |\n")

        f.write("\n")

        # Domain-specific confidence
        f.write("## Domain-Specific Confidence (Mean pLDDT)\n\n")
        f.write("| Domain | Model 0 | Model 1 | Model 2 | Model 3 | Model 4 | Status |\n")
        f.write("|--------|---------|---------|---------|---------|---------|--------|\n")

        for domain_name in DOMAINS.keys():
            f.write(f"| {domain_name:20s} | ")
            scores = []
            for result in results:
                if domain_name in result['domains']:
                    score = result['domains'][domain_name]['mean']
                    scores.append(score)
                    f.write(f"{score:7.2f} | ")
                else:
                    f.write("    N/A | ")

            # Status
            if scores:
                avg_score = np.mean(scores)
                if avg_score >= 90:
                    status = "✅ Excellent"
                elif avg_score >= 70:
                    status = "✅ Good"
                elif avg_score >= 50:
                    status = "⚠️ Low"
                else:
                    status = "❌ Very Low"
                f.write(f" {status} |\n")
            else:
                f.write(" N/A |\n")

        f.write("\n")

        # Detailed domain statistics for best model (Model 0)
        f.write("## Detailed Domain Statistics (Model 0 - Best Ranked)\n\n")
        f.write("| Domain | Mean | Median | Min | Max | Std | Residues |\n")
        f.write("|--------|------|--------|-----|-----|-----|----------|\n")

        for domain_name, stats in results[0]['domains'].items():
            f.write(f"| {domain_name:20s} | "
                   f"{stats['mean']:6.2f} | "
                   f"{stats['median']:6.2f} | "
                   f"{stats['min']:6.2f} | "
                   f"{stats['max']:6.2f} | "
                   f"{stats['std']:6.2f} | "
                   f"{stats['num_residues']:8d} |\n")

        f.write("\n")

        # Recommendations
        f.write("## Quality Assessment\n\n")

        best_model = results[0]
        hect_plddt = best_model['domains']['HECT']['mean']
        ww_plddt = (best_model['domains']['WW1']['mean'] + best_model['domains']['WW2']['mean']) / 2
        c2_plddt = best_model['domains']['C2']['mean']

        f.write(f"**HECT Domain (aa 1268-1571):** {hect_plddt:.2f}\n")
        if hect_plddt >= 70:
            f.write("- ✅ PASS: Sufficient confidence for FoldX/Rosetta analysis\n")
        else:
            f.write("- ⚠️ WARNING: Low confidence, results may be unreliable\n")
        f.write("\n")

        f.write(f"**WW Domains (aa 809-838, 987-1016):** {ww_plddt:.2f}\n")
        if ww_plddt >= 60:
            f.write("- ✅ PASS: Acceptable confidence for small domains\n")
        else:
            f.write("- ⚠️ WARNING: Low confidence, WW domains may be poorly modeled\n")
        f.write("\n")

        f.write(f"**C2 Domain (aa 192-281):** {c2_plddt:.2f}\n")
        if c2_plddt >= 70:
            f.write("- ✅ PASS: Good confidence\n")
        else:
            f.write("- ⚠️ WARNING: Low confidence\n")
        f.write("\n")

        f.write("**Overall Assessment:**\n")
        if hect_plddt >= 70 and best_model['overall']['has_clash'] == 0:
            f.write("- ✅ **READY FOR PHASE 3**: Structures suitable for mutagenesis and ΔΔG analysis\n")
        elif hect_plddt >= 60:
            f.write("- ⚠️ **PROCEED WITH CAUTION**: HECT domain confidence acceptable but not ideal\n")
        else:
            f.write("- ❌ **NOT RECOMMENDED**: Low confidence may affect ΔΔG accuracy\n")

        f.write("\n")
        f.write("## Next Steps\n\n")
        f.write("1. Select 3 diverse seeds from 5 models (RMSD analysis)\n")
        f.write("2. Verify no missing residues in HECT domain\n")
        f.write("3. Proceed to Phase 3: In silico mutagenesis\n")

def main():
    # Paths
    base_dir = '/Users/ykshim2025/Desktop/Code2025/HECW2/fold_hecw2_2025_10_10_18_12'
    output_dir = '/Users/ykshim2025/Desktop/Code2025/HECW2/structures/alphafold_output'
    quality_dir = '/Users/ykshim2025/Desktop/Code2025/HECW2/structures/quality'

    # Analyze all 5 models
    print("Analyzing AlphaFold3 output for HECW2...")
    print("="*70)

    results = []
    for model_num in range(5):
        print(f"\nProcessing Model {model_num}...")
        result = analyze_model(model_num, base_dir, output_dir)
        results.append(result)
        print(f"  Mean pLDDT: {result['overall']['mean_plddt']:.2f}")
        print(f"  HECT domain pLDDT: {result['domains']['HECT']['mean']:.2f}")

    # Generate report
    report_file = os.path.join(quality_dir, 'alphafold3_quality_report.txt')
    generate_quality_report(results, report_file)
    print(f"\n✅ Quality report saved: {report_file}")

    # Save pLDDT scores to CSV
    for result in results:
        csv_file = os.path.join(quality_dir, f"plddt_model{result['model']}.csv")
        with open(csv_file, 'w') as f:
            f.write("Residue,pLDDT\n")
            for res_id, plddt in result['plddt_per_residue']:
                f.write(f"{res_id},{plddt:.2f}\n")
        print(f"✅ pLDDT scores saved: {csv_file}")

    print("\n" + "="*70)
    print("Analysis complete!")
    print(f"PDB files: {output_dir}")
    print(f"Quality report: {report_file}")

if __name__ == '__main__':
    main()
