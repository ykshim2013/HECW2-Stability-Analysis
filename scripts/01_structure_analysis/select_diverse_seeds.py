#!/usr/bin/env python3
"""
Select 3 diverse seeds from 5 AlphaFold3 models
- Calculate pairwise RMSD
- Select models with maximum diversity
- Verify no missing residues in HECT domain
"""

import os
import numpy as np
from Bio.PDB import PDBParser, Superimposer, PDBIO
import shutil

def calculate_rmsd(structure1, structure2, domain_range=None):
    """Calculate RMSD between two structures"""
    # Get CA atoms
    atoms1 = []
    atoms2 = []

    for model in structure1:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':  # Standard residue
                    if domain_range:
                        res_id = residue.id[1]
                        if not (domain_range[0] <= res_id <= domain_range[1]):
                            continue
                    ca_atoms = [atom for atom in residue if atom.name == 'CA']
                    if ca_atoms:
                        atoms1.append(ca_atoms[0])

    for model in structure2:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':
                    if domain_range:
                        res_id = residue.id[1]
                        if not (domain_range[0] <= res_id <= domain_range[1]):
                            continue
                    ca_atoms = [atom for atom in residue if atom.name == 'CA']
                    if ca_atoms:
                        atoms2.append(ca_atoms[0])

    # Ensure same number of atoms
    if len(atoms1) != len(atoms2):
        print(f"Warning: Atom count mismatch ({len(atoms1)} vs {len(atoms2)})")
        min_len = min(len(atoms1), len(atoms2))
        atoms1 = atoms1[:min_len]
        atoms2 = atoms2[:min_len]

    # Calculate RMSD
    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)
    return super_imposer.rms

def verify_hect_completeness(structure, hect_range=(1268, 1571)):
    """Verify no missing residues in HECT domain"""
    residue_ids = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':
                    res_id = residue.id[1]
                    if hect_range[0] <= res_id <= hect_range[1]:
                        residue_ids.append(res_id)

    residue_ids = sorted(set(residue_ids))
    expected = list(range(hect_range[0], hect_range[1] + 1))
    missing = set(expected) - set(residue_ids)

    return len(missing) == 0, missing

def select_diverse_models(pdb_dir, n_models=3):
    """Select n most diverse models based on RMSD"""
    parser = PDBParser(QUIET=True)

    # Load all models
    models = []
    for i in range(5):
        pdb_file = os.path.join(pdb_dir, f'HECW2_model{i}.pdb')
        structure = parser.get_structure(f'model{i}', pdb_file)
        models.append((i, structure))

    # Calculate pairwise RMSD matrix (HECT domain only)
    n = len(models)
    rmsd_matrix = np.zeros((n, n))

    print("Calculating pairwise RMSD (HECT domain aa 1268-1571)...")
    print("="*70)

    for i in range(n):
        for j in range(i+1, n):
            rmsd = calculate_rmsd(models[i][1], models[j][1], domain_range=(1268, 1571))
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd
            print(f"Model {i} vs Model {j}: RMSD = {rmsd:.2f} Å")

    print("\n" + "="*70)
    print("RMSD Matrix (HECT domain):\n")
    print("       ", end="")
    for i in range(n):
        print(f"Model{i:1d}  ", end="")
    print()

    for i in range(n):
        print(f"Model{i} ", end="")
        for j in range(n):
            if i == j:
                print("  -    ", end="")
            else:
                print(f"{rmsd_matrix[i, j]:5.2f}  ", end="")
        print()

    # Select most diverse set using greedy algorithm
    # Start with model 0 (best ranked by AlphaFold3)
    selected = [0]

    # Iteratively add model with maximum minimum RMSD to selected set
    while len(selected) < n_models:
        max_min_rmsd = -1
        best_candidate = None

        for candidate in range(n):
            if candidate in selected:
                continue

            # Calculate minimum RMSD to already selected models
            min_rmsd = min(rmsd_matrix[candidate, s] for s in selected)

            if min_rmsd > max_min_rmsd:
                max_min_rmsd = min_rmsd
                best_candidate = candidate

        selected.append(best_candidate)

    print("\n" + "="*70)
    print(f"Selected {n_models} diverse seeds:")
    for i, model_idx in enumerate(selected, 1):
        print(f"  Seed {i}: Model {model_idx}")

    # Calculate diversity metrics
    selected_rmsds = []
    for i in range(len(selected)):
        for j in range(i+1, len(selected)):
            selected_rmsds.append(rmsd_matrix[selected[i], selected[j]])

    print(f"\nDiversity metrics:")
    print(f"  Mean pairwise RMSD: {np.mean(selected_rmsds):.2f} Å")
    print(f"  Min pairwise RMSD:  {np.min(selected_rmsds):.2f} Å")
    print(f"  Max pairwise RMSD:  {np.max(selected_rmsds):.2f} Å")

    return selected, models

def main():
    pdb_dir = '/Users/ykshim2025/Desktop/Code2025/HECW2/structures/alphafold_output'
    seed_dir = '/Users/ykshim2025/Desktop/Code2025/HECW2/structures/seeds'

    print("HECW2 Seed Selection for Phase 3")
    print("="*70)
    print()

    # Select diverse models
    selected_indices, models = select_diverse_models(pdb_dir, n_models=3)

    # Verify HECT domain completeness
    print("\n" + "="*70)
    print("Verifying HECT domain completeness (aa 1268-1571)...")
    print("="*70)

    all_complete = True
    for seed_num, model_idx in enumerate(selected_indices, 1):
        is_complete, missing = verify_hect_completeness(models[model_idx][1])

        if is_complete:
            print(f"✅ Seed {seed_num} (Model {model_idx}): HECT domain complete (304 residues)")
        else:
            print(f"❌ Seed {seed_num} (Model {model_idx}): Missing {len(missing)} residues: {sorted(missing)}")
            all_complete = False

    if not all_complete:
        print("\n⚠️ WARNING: Some models have missing residues in HECT domain!")
        print("This may affect ΔΔG calculations for variants in missing regions.")
        return

    # Copy selected models to seeds directory
    print("\n" + "="*70)
    print(f"Copying selected models to {seed_dir}...")
    print("="*70)

    os.makedirs(seed_dir, exist_ok=True)

    for seed_num, model_idx in enumerate(selected_indices, 1):
        src = os.path.join(pdb_dir, f'HECW2_model{model_idx}.pdb')
        dst = os.path.join(seed_dir, f'seed{seed_num}.pdb')
        shutil.copy2(src, dst)
        print(f"✅ Copied Model {model_idx} → seed{seed_num}.pdb")

    print("\n" + "="*70)
    print("✅ SEED SELECTION COMPLETE!")
    print("="*70)
    print(f"\nSelected seeds saved to: {seed_dir}")
    print("Ready for Phase 3: In silico mutagenesis")
    print("\nNext steps:")
    print("  1. Review seed diversity metrics")
    print("  2. Verify PDB files in structures/seeds/")
    print("  3. Proceed to Phase 3: Generate 31 × 3 = 93 mutant structures")

if __name__ == '__main__':
    main()
