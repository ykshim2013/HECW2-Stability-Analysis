# Supplementary Methods

## HECW2 Variant Stability Analysis: Detailed Protocols for Reproducibility

---

## Table of Contents

1. [Software and Dependencies](#1-software-and-dependencies)
2. [AlphaFold3 Structure Prediction](#2-alphafold3-structure-prediction)
3. [Seed Selection Algorithm](#3-seed-selection-algorithm)
4. [FoldX Analysis Protocol](#4-foldx-analysis-protocol)
5. [Rosetta cartesian_ddG Protocol](#5-rosetta-cartesian_ddg-protocol)
6. [Result Parsing and Statistical Analysis](#6-result-parsing-and-statistical-analysis)
7. [Data Files and Repository Structure](#7-data-files-and-repository-structure)

---

## 1. Software and Dependencies

### 1.1 Software Versions

| Software | Version | Purpose |
|----------|---------|---------|
| AlphaFold3 | Server (October 2024) | Structure prediction |
| FoldX | 5.1 | Stability prediction (BuildModel + Stability) |
| Rosetta | 2023.49 | Stability prediction (cartesian_ddG) |
| Python | 3.10+ | Analysis scripts |
| BioPython | 1.83 | PDB parsing and manipulation |
| NumPy | 1.24+ | Numerical computation |
| Pandas | 2.0+ | Data manipulation |
| SciPy | 1.11+ | Statistical analysis |
| Matplotlib | 3.7+ | Visualization |

### 1.2 Python Dependencies

```bash
# Install required packages
pip install biopython numpy pandas scipy matplotlib seaborn scikit-learn
```

### 1.3 System Requirements

- **Operating System**: macOS 12+ or Linux (Ubuntu 20.04+)
- **RAM**: Minimum 16 GB (32 GB recommended for Rosetta)
- **CPU**: Multi-core processor (8+ cores recommended for parallel processing)
- **Disk Space**: ~50 GB for complete analysis including all intermediate files

### 1.4 Variant Selection Criteria

Variants were selected to ensure representation across all HECW2 functional domains while prioritizing those with phenotype descriptions available in ClinVar. Synonymous, intronic, UTR, and frameshift variants were excluded as these variant types cannot be modeled by ΔΔG analysis. A total of 96 missense variants were analyzed: 34 pathogenic/likely pathogenic, 22 benign/likely benign, and 40 variants of uncertain significance (VUS).

---

## 2. AlphaFold3 Structure Prediction

### 2.1 Input Sequence

Full-length HECW2 sequence (UniProt Q96JN2, 1,572 amino acids) was submitted to the AlphaFold3 Server. Complete full-length AlphaFold3-predicted structures were used as input for both FoldX and Rosetta methods, rather than isolated domains, to preserve the native structural context and inter-domain interactions relevant to variant impact.

**Input file**: `alphafold_sequence.fasta`
```
>HECW2_Human|UniProt:Q96JN2
MAQVQSQSEG... [full 1,572 aa sequence]
```

### 2.2 Server Parameters

- **Server**: AlphaFold3 Server (https://alphafoldserver.com)
- **Date**: October 2024
- **Template Constraints**: None (de novo prediction)
- **Number of Models**: 5 (default)
- **Relaxation**: Standard AlphaFold3 relaxation

### 2.3 Output Files

AlphaFold3 generates the following output for each prediction:

```
fold_hecw2_2025_10_10_18_12/
├── fold_hecw2_model_0.cif          # Model 0 (rank 1)
├── fold_hecw2_model_1.cif          # Model 1 (rank 2)
├── fold_hecw2_model_2.cif          # Model 2 (rank 3)
├── fold_hecw2_model_3.cif          # Model 3 (rank 4)
├── fold_hecw2_model_4.cif          # Model 4 (rank 5)
├── fold_hecw2_full_data_0.json     # Confidence data for model 0
├── fold_hecw2_full_data_1.json     # Confidence data for model 1
├── fold_hecw2_full_data_2.json     # Confidence data for model 2
├── fold_hecw2_full_data_3.json     # Confidence data for model 3
├── fold_hecw2_full_data_4.json     # Confidence data for model 4
└── fold_hecw2_job_request.json     # Job parameters
```

### 2.4 CIF to PDB Conversion

```bash
# Convert mmCIF to PDB format using BioPython
python3 scripts/convert_cif_to_pdb.py
```

### 2.5 pLDDT Score Extraction

Per-residue predicted Local Distance Difference Test (pLDDT) scores were extracted directly from AlphaFold3 confidence output files to assess structural quality at variant positions. For each variant position, pLDDT values were averaged across the three structural seeds and reported with standard deviation. Domain-specific quality assessment was performed with pLDDT >70 established as the reliability threshold for energy calculations.

Per-residue pLDDT scores were extracted from the JSON confidence files:

```python
import json

def extract_plddt(json_file):
    """Extract per-residue pLDDT scores from AlphaFold3 output"""
    with open(json_file, 'r') as f:
        data = json.load(f)

    # pLDDT is stored in the 'atom_plddts' field
    plddt_scores = data['atom_plddts']

    # Average per residue (multiple atoms per residue)
    # Implementation depends on AlphaFold3 output format
    return plddt_scores
```

---

## 3. Seed Selection Algorithm

### 3.1 Rationale

Three conformational seeds were selected from 5 AlphaFold3 models to:
1. Capture structural ensemble diversity
2. Reduce prediction variance
3. Follow established protocols (FBXO11 precedent)

### 3.2 Selection Algorithm

**Script**: `scripts/select_diverse_seeds.py`

```python
def select_diverse_models(pdb_dir, n_models=3):
    """
    Select n most diverse models based on RMSD using greedy algorithm.

    Algorithm:
    1. Start with model 0 (best AlphaFold3 ranking)
    2. Iteratively add model with maximum minimum RMSD to selected set
    3. Repeat until n models selected
    """
    # Calculate pairwise RMSD matrix (HECT domain only, aa 1268-1571)
    rmsd_matrix = calculate_pairwise_rmsd(models, domain_range=(1268, 1571))

    # Greedy selection
    selected = [0]  # Start with best-ranked model

    while len(selected) < n_models:
        max_min_rmsd = -1
        best_candidate = None

        for candidate in range(len(models)):
            if candidate in selected:
                continue

            min_rmsd = min(rmsd_matrix[candidate, s] for s in selected)

            if min_rmsd > max_min_rmsd:
                max_min_rmsd = min_rmsd
                best_candidate = candidate

        selected.append(best_candidate)

    return selected
```

### 3.3 RMSD Calculation

RMSD was calculated using Cα atoms only, restricted to the HECT domain (residues 1268–1571):

```python
from Bio.PDB import Superimposer

def calculate_rmsd(structure1, structure2, domain_range=(1268, 1571)):
    """Calculate Cα RMSD between two structures for specified domain"""
    # Extract Cα atoms within domain range
    atoms1 = [atom for residue in structure1
              if domain_range[0] <= residue.id[1] <= domain_range[1]
              for atom in residue if atom.name == 'CA']

    atoms2 = [atom for residue in structure2
              if domain_range[0] <= residue.id[1] <= domain_range[1]
              for atom in residue if atom.name == 'CA']

    # Superimpose and calculate RMSD
    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)
    return super_imposer.rms
```

### 3.4 Selected Seeds

| Seed | AlphaFold3 Model | Mean RMSD to Others |
|------|------------------|---------------------|
| Seed 1 | Model 0 (rank 1) | 0.71 Å |
| Seed 2 | Model 1 (rank 2) | 0.68 Å |
| Seed 3 | Model 4 (rank 5) | 0.79 Å |

**Mean pairwise RMSD**: 0.73 Å (indicating converged predictions)

---

## 4. FoldX Analysis Protocol

### 4.1 Overview

FoldX 5.1 was used to calculate ΔΔG values using the BuildModel command with 10 technical replicates per variant per seed.

### 4.2 Structure Preparation (RepairPDB)

Before mutagenesis, each seed structure was optimized using FoldX RepairPDB:

```bash
#!/bin/bash
# Repair wildtype structures to optimize hydrogen bonding and remove clashes

FOLDX="/path/to/foldx"

for seed in seed1 seed2 seed3; do
    $FOLDX --command=RepairPDB \
           --pdb=${seed}.pdb \
           --output-dir=repaired_wt/
done
```

**Purpose of RepairPDB**:
- Optimize hydrogen bonding networks
- Remove steric clashes
- Energy minimize local regions
- Ensure consistent starting point for ΔΔG calculations

### 4.3 Mutation File Format

FoldX uses a specific format for mutation files (`individual_list.txt`):

```
# Format: [WT_AA][Chain][Position][Mut_AA];
RA15Q;
IA138V;
EA39A;
TA69A;
EA744K;
DA1134G;
PA1151S;
...
```

**Note**: Chain identifier 'A' is included for single-chain structures.

### 4.4 BuildModel Execution

**Script**: `scripts/run_foldx_pipeline.sh`

```bash
#!/bin/bash
# FoldX BuildModel for all variants on all seeds

FOLDX="/path/to/foldx"
N_RUNS=10  # Technical replicates

for seed in seed1 seed2 seed3; do
    # BuildModel generates mutant structures and calculates ΔΔG
    $FOLDX --command=BuildModel \
           --pdb=wt_${seed}.pdb \
           --mutant-file=individual_list.txt \
           --numberOfRuns=$N_RUNS \
           --output-dir=./${seed}_foldx/
done
```

### 4.5 FoldX Output Files

```
seed1_foldx/
├── Dif_wt_seed1.fxout           # ΔΔG values (primary output)
├── Average_wt_seed1.fxout       # Average energies
├── Raw_wt_seed1.fxout           # Detailed energy terms
├── wt_seed1_1_0.pdb             # Mutant structure (variant 1, run 0)
├── wt_seed1_1_1.pdb             # Mutant structure (variant 1, run 1)
├── ...
└── buildmodel.log               # Execution log
```

### 4.6 ΔΔG Calculation and Threshold Rationale

FoldX calculates ΔΔG as:

```
ΔΔG = ΔG_mutant - ΔG_wildtype
```

Where ΔG is the total energy from the FoldX force field.

**Output format** (`Dif_wt_seed1.fxout`):
```
Pdb	total energy	Backbone Hbond	Sidechain Hbond	...
wt_seed1_1_0.pdb	-0.102171	0.523728	-0.456123	...
wt_seed1_1_1.pdb	-0.105994	-0.603671	0.234567	...
```

**ΔΔG Threshold Rationale**: ΔΔG thresholds of 2 kcal/mol and 5 kcal/mol were employed based on established literature. The 2 kcal/mol threshold discriminates neutral from perturbing substitutions and has been widely adopted in structure-based variant analysis. The 5 kcal/mol threshold identifies severely destabilizing variants, as previous benchmarks demonstrated that thresholds in this range yield high positive predictive value (>95%) for identifying functionally damaging variants. These thresholds also exceed typical computational uncertainty (±1–2 kcal/mol).

### 4.7 Runtime Estimates

| Mode | Variants | Seeds | Runs | Total Calculations | Time |
|------|----------|-------|------|-------------------|------|
| Test | 3 | 1 | 10 | 30 | ~10 min |
| Full | 96 | 3 | 10 | 2,880 | ~15 hours |

---

## 5. Rosetta cartesian_ddG Protocol

### 5.1 Overview

Rosetta cartesian_ddG was used as an orthogonal validation method with the ref2015_cart energy function and backrub ensemble sampling.

### 5.2 Structure Relaxation (FastRelax)

Before ΔΔG calculations, structures were relaxed in Rosetta:

```bash
#!/bin/bash
# Relax structures to Rosetta energy function

ROSETTA_BIN="/path/to/rosetta/source/bin"
ROSETTA_DB="/path/to/rosetta/database"

for seed in seed1 seed2 seed3; do
    ${ROSETTA_BIN}/relax.macosclangrelease \
        -database ${ROSETTA_DB} \
        -s ${seed}.pdb \
        -relax:constrain_relax_to_start_coords \
        -relax:coord_constrain_sidechains \
        -relax:ramp_constraints false \
        -ex1 -ex2 \
        -use_input_sc \
        -flip_HNQ \
        -no_optH false \
        -out:prefix relaxed_ \
        -nstruct 10
done
```

**Parameters explained**:
- `-relax:constrain_relax_to_start_coords`: Maintain overall structure
- `-relax:coord_constrain_sidechains`: Constrain sidechain positions
- `-ex1 -ex2`: Expanded rotamer sampling
- `-nstruct 10`: Generate 10 relaxed background structures

### 5.3 Mutation File Format

Rosetta uses a different format for mutations:

```
# Format: total N_mutations
# N_mutations
# [WT_AA] [Position] [Chain] PIKAA [Mut_AA]

total 1
1
R 15 A PIKAA Q
```

### 5.4 cartesian_ddG Execution

**Script**: `scripts/run_rosetta_cartesian_ddg.sh`

```bash
#!/bin/bash
# Rosetta cartesian_ddG for all variants

ROSETTA_BIN="/path/to/rosetta/source/bin/cartesian_ddg.macosclangrelease"
ROSETTA_DB="/path/to/rosetta/database"
NUM_ITERATIONS=17

for variant_mut in mutation_files/*.mut; do
    for bg_pdb in relaxed_backgrounds/*.pdb; do
        ${ROSETTA_BIN} \
            -database ${ROSETTA_DB} \
            -s ${bg_pdb} \
            -ddg:mut_file ${variant_mut} \
            -ddg:iterations ${NUM_ITERATIONS} \
            -ddg:cartesian \
            -score:weights ref2015_cart \
            -fa_max_dis 9.0 \
            -ddg:bbnbr 1 \
            -ddg:dump_pdbs false \
            -ddg:suppress_checkpointing true \
            -ddg:mean false \
            -ddg:min true \
            -ddg:sc_min_only false \
            -ddg:opt_radius 6.0 \
            -mute all
    done
done
```

### 5.5 Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `-ddg:cartesian` | true | Use Cartesian space optimization |
| `-score:weights` | ref2015_cart | Cartesian-optimized energy function |
| `-ddg:iterations` | 17 | Backrub sampling iterations |
| `-ddg:opt_radius` | 6.0 | Optimization radius around mutation site (Å) |
| `-ddg:bbnbr` | 1 | Allow backbone flexibility (1 neighbor) |
| `-fa_max_dis` | 9.0 | Maximum distance for full-atom interactions |

### 5.6 Rosetta Output Files

```
results/p.ARG15GLN_seed1/
├── mutations.ddg              # ΔΔG predictions
├── p.ARG15GLN_seed1.log       # Execution log
└── ddg_predictions.out        # Formatted output
```

**Output format** (`mutations.ddg`):
```
COMPLEX:   Round1: WT: 15267.858  fa_atr: -8234.123 ...
COMPLEX:   Round1: MUT_15GLN: 15272.345  fa_atr: -8229.456 ...
COMPLEX:   Round2: WT: 15268.123  fa_atr: -8235.789 ...
COMPLEX:   Round2: MUT_15GLN: 15273.567  fa_atr: -8230.012 ...
...
```

### 5.7 ΔΔG Calculation

```python
def calculate_rosetta_ddg(wt_energies, mut_energies):
    """
    Calculate ΔΔG = <E_mut> - <E_wt>

    Error propagation:
    σ_ΔΔG = sqrt(σ_mut² + σ_wt²)
    """
    mean_wt = np.mean(wt_energies)
    mean_mut = np.mean(mut_energies)

    ddg = mean_mut - mean_wt

    sd_wt = np.std(wt_energies, ddof=1)
    sd_mut = np.std(mut_energies, ddof=1)
    sd_ddg = np.sqrt(sd_wt**2 + sd_mut**2)

    return ddg, sd_ddg
```

### 5.8 Runtime Estimates

| Configuration | Time per Variant | Total (96 variants × 3 seeds) |
|---------------|------------------|-------------------------------|
| Sequential | ~30 min | ~140 hours |
| Parallel (10 cores) | ~3 min | ~14 hours |

---

## 6. Result Parsing and Statistical Analysis

### 6.1 FoldX Result Parsing

**Script**: `scripts/parse_foldx_results.py`

```python
def parse_foldx_dif_file(filepath):
    """
    Parse FoldX Dif_*.fxout file

    Returns DataFrame with columns:
    - pdb_file, ddg_kcal_mol, seed, variant_id, run_id
    """
    # Read file, skip header lines
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Find header line (starts with "Pdb")
    header_idx = next(i for i, line in enumerate(lines)
                      if line.startswith('Pdb'))

    # Parse data lines
    data = []
    for line in lines[header_idx + 1:]:
        parts = line.strip().split('\t')
        pdb_file = parts[0]
        ddg = float(parts[1])  # Total energy (ΔΔG)

        # Extract variant_id and run_id from filename
        data.append({
            'pdb_file': pdb_file,
            'ddg_kcal_mol': ddg,
            ...
        })

    return pd.DataFrame(data)
```

### 6.2 Rosetta Result Parsing

**Script**: `scripts/parse_rosetta_ddg.py`

```python
def parse_rosetta_ddg_file(filepath):
    """
    Parse Rosetta mutations.ddg file

    Returns lists of WT and mutant energies per round
    """
    wt_energies = []
    mut_energies = []

    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith('COMPLEX:'):
                continue

            # Extract WT energy
            match_wt = re.search(r'WT:\s+([\d.]+)', line)
            if match_wt:
                wt_energies.append(float(match_wt.group(1)))

            # Extract mutant energy
            match_mut = re.search(r'MUT_\w+:\s+([\d.]+)', line)
            if match_mut:
                mut_energies.append(float(match_mut.group(1)))

    return wt_energies, mut_energies
```

### 6.3 Statistical Aggregation

**FoldX**: Median ± IQR (robust to outliers)
```python
def aggregate_foldx(df):
    """Aggregate FoldX results using median ± IQR"""
    grouped = df.groupby('variant')['ddg_kcal_mol']

    summary = pd.DataFrame({
        'median_ddg': grouped.median(),
        'iqr': grouped.apply(lambda x: x.quantile(0.75) - x.quantile(0.25)),
        'n_measurements': grouped.count()
    })

    return summary
```

**Rosetta**: Mean ± SD (standard for backrub sampling)
```python
def aggregate_rosetta(df):
    """Aggregate Rosetta results using mean ± SD"""
    grouped = df.groupby('variant')['ddg_kcal_mol']

    summary = pd.DataFrame({
        'mean_ddg': grouped.mean(),
        'sd': grouped.std(ddof=1),
        'n_seeds': grouped.count()
    })

    return summary
```

### 6.4 Cross-Seed Aggregation

Final variant-level ΔΔG values were calculated by averaging across three seeds:

```python
def aggregate_across_seeds(df):
    """
    For each variant:
    1. Calculate seed-level statistics (mean of technical replicates)
    2. Average across 3 seeds
    3. Report mean ± SD across seeds
    """
    # Step 1: Seed-level means
    seed_means = df.groupby(['variant', 'seed'])['ddg_kcal_mol'].mean()

    # Step 2: Variant-level statistics
    variant_stats = seed_means.groupby('variant').agg(['mean', 'std', 'count'])

    return variant_stats
```

### 6.5 Statistical Tests

**Discrimination analysis**:
```python
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve

# Mann-Whitney U test (non-parametric)
stat, pvalue = stats.mannwhitneyu(
    pathogenic_ddg, benign_ddg,
    alternative='greater'
)

# Cohen's d effect size
cohens_d = (pathogenic_ddg.mean() - benign_ddg.mean()) / pooled_std

# ROC analysis
auc = roc_auc_score(y_true, ddg_values)
fpr, tpr, thresholds = roc_curve(y_true, ddg_values)

# 95% CI for AUC (DeLong method)
auc_ci = delong_ci(y_true, ddg_values)
```

---

## 7. Data Files and Repository Structure

### 7.1 GitHub Repository Structure

```
HECW2-stability-analysis/
├── README.md
├── LICENSE
├── requirements.txt
│
├── data/
│   ├── variants/
│   │   ├── variant_list.csv           # All 96 variants with classifications
│   │   ├── clinvar_annotations.csv    # ClinVar data
│   │   └── phenotype_data.csv         # Clinical phenotypes
│   ├── structures/
│   │   ├── alphafold_models/          # AlphaFold3 output (CIF)
│   │   └── seeds/                     # Selected seed PDBs
│   └── results/
│       ├── foldx_ddg_summary.csv      # FoldX results
│       ├── rosetta_ddg_summary.csv    # Rosetta results
│       └── combined_ddg_analysis.csv  # Dual-method comparison
│
├── scripts/
│   ├── structure_analysis/
│   │   ├── analyze_alphafold_output.py
│   │   └── select_diverse_seeds.py
│   ├── foldx/
│   │   ├── run_foldx_pipeline.sh
│   │   └── parse_foldx_results.py
│   ├── rosetta/
│   │   ├── run_rosetta_cartesian_ddg.sh
│   │   └── parse_rosetta_ddg.py
│   └── analysis/
│       ├── combine_rosetta_foldx.py
│       └── statistical_analysis.py
│
├── manuscript/
│   ├── figures/
│   └── tables/
│
└── docs/
    ├── Supplementary_Methods.md
    └── Supplementary_Tables.xlsx
```

### 7.2 Data File Formats

**variant_list.csv**:
```csv
variant,clinical_classification,domain,position,wt_aa,mut_aa,clinvar_id
p.Arg15Gln,VUS,N-terminal,15,R,Q,VCV001234567
p.Ile138Val,Likely_Pathogenic,Inter1,138,I,V,VCV001234568
...
```

**foldx_ddg_summary.csv**:
```csv
variant,median_ddg_kcal_mol,iqr,mean_ddg_kcal_mol,sd,n_measurements,classification
p.Gly1452Val,17.68,2.34,18.12,3.45,30,Pathogenic
p.Ser1504Tyr,12.45,1.89,12.78,2.67,30,Pathogenic
...
```

**rosetta_ddg_summary.csv**:
```csv
variant,mean_ddg_kcal_mol,sd_across_seeds,mean_technical_sd,n_seeds,classification
p.Gly1452Val,15.23,2.12,1.45,3,Pathogenic
p.Ser1504Tyr,10.89,1.78,1.23,3,Pathogenic
...
```

### 7.3 Reproducibility Checklist

- [ ] Clone repository: `git clone https://github.com/ykshim2013/HECW2-Stability-Analysis.git`
- [ ] Install dependencies: `pip install -r requirements.txt`
- [ ] Download AlphaFold3 structures (if not included)
- [ ] Install FoldX 5.1 (requires academic license)
- [ ] Install Rosetta 2023.49 (requires academic license)
- [ ] Run FoldX pipeline: `./scripts/foldx/run_foldx_pipeline.sh full`
- [ ] Run Rosetta pipeline: `./scripts/rosetta/run_rosetta_cartesian_ddg.sh`
- [ ] Parse results: `python3 scripts/analysis/combine_rosetta_foldx.py`

---

## References

1. Delgado J, et al. FoldX 5.0: working with RNA, small molecules and a new graphical interface. *Bioinformatics*. 2019;35(20):4168-4169.

2. Park H, et al. Simultaneous optimization of biomolecular energy functions on features from small molecules and macromolecules. *J Chem Theory Comput*. 2016;12(12):6201-6212.

3. Abramson J, et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature*. 2024;630:493-500.

4. Guerois R, et al. Predicting changes in the stability of proteins and protein complexes: a study of more than 1000 mutations. *J Mol Biol*. 2002;320(2):369-387.

5. Kellogg EH, et al. Role of conformational sampling in computing mutation-induced changes in protein structure and stability. *Proteins*. 2011;79(3):830-838.

---

*Last updated: December 2024*

*Corresponding author:*
Youngkyu Shim, M.D.
Department of Pediatrics, Korea University Ansan Hospital,
Korea University College of Medicine
