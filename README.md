# HECW2 Variant Stability Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Computational pipeline for dual-method ΔΔG stability predictions (FoldX + Rosetta) of HECW2 missense variants for variant interpretation in Neurodevelopmental Disorder with Hypotonia, Seizures, and Absent Language (NDHSAL).

---

## Overview

This repository contains the computational pipeline and analysis scripts for reproducible stability analysis:

- **AlphaFold3** structure prediction and quality assessment
- **FoldX 5.1** stability predictions (BuildModel)
- **Rosetta 2023.49** stability predictions (cartesian_ddG)

**Associated Publication:**
> Shim YK. HECW2 Pathogenic Variants Exhibit Mechanistic Heterogeneity Revealed by Dual-Method Computational Stability Analysis. *Frontiers in Neurology*. 2025. (in submission)

---

## Repository Structure

```
HECW2-stability-analysis/
├── README.md
├── LICENSE
├── requirements.txt
├── Supplementary_Methods.md
│
├── data/
│   └── structures/
│       └── seeds/                    # 3 AlphaFold3 seed structures (PDB)
│           ├── seed1.pdb
│           ├── seed2.pdb
│           ├── seed3.pdb
│           ├── seed1_hect_domain.pdb
│           └── seed1_hect_core.pdb
│
├── scripts/
│   ├── 01_structure_analysis/
│   │   ├── analyze_alphafold_output.py   # pLDDT extraction and QC
│   │   └── select_diverse_seeds.py       # RMSD-based seed selection
│   │
│   ├── 02_foldx/
│   │   ├── run_foldx_pipeline.sh         # Main FoldX workflow
│   │   └── parse_foldx_results.py        # Parse Dif_*.fxout files
│   │
│   ├── 03_rosetta/
│   │   ├── run_rosetta_cartesian_ddg.sh  # Main Rosetta workflow
│   │   └── parse_rosetta_ddg.py          # Parse mutations.ddg files
│   │
│   └── 04_analysis/
│       └── combine_rosetta_foldx.py      # Merge dual-method results
│
└── Supplementary_Methods.md              # Detailed protocols
```

---

## Quick Start

### 1. Clone Repository

```bash
git clone https://github.com/ykshim2013/HECW2-Stability-Analysis.git
cd HECW2-Stability-Analysis
```

### 2. Install Dependencies

```bash
# Create virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate

# Install Python packages
pip install -r requirements.txt
```

### 3. External Software Requirements

| Software | Version | License | Download |
|----------|---------|---------|----------|
| FoldX | 5.1 | Academic | [FoldX](https://foldxsuite.crg.eu/) |
| Rosetta | 2023.49 | Academic | [RosettaCommons](https://www.rosettacommons.org/) |

### 4. Run Analysis Pipeline

```bash
# Step 1: Seed selection (if starting from AlphaFold3 output)
python3 scripts/01_structure_analysis/select_diverse_seeds.py

# Step 2: FoldX analysis
./scripts/02_foldx/run_foldx_pipeline.sh full

# Step 3: Rosetta analysis
./scripts/03_rosetta/run_rosetta_cartesian_ddg.sh

# Step 4: Parse and combine results
python3 scripts/04_analysis/combine_rosetta_foldx.py
```

---

## Methods Summary

### AlphaFold3 Structure Prediction
- Full-length HECW2 (1,572 aa) predicted using AlphaFold3 Server (October 2024)
- 5 models generated, 3 diverse seeds selected based on HECT domain RMSD
- Mean pLDDT for HECT domain: 88.5 (high confidence)

### FoldX Protocol
1. **RepairPDB**: Optimize H-bonding and remove clashes
2. **BuildModel**: Generate mutant structures (10 replicates × 3 seeds)
3. **Statistics**: Mean ± SD across seeds

### Rosetta Protocol
1. **FastRelax**: Energy minimize in Rosetta force field
2. **cartesian_ddG**: 17 iterations with ref2015_cart energy function
3. **Statistics**: Mean ± SD across seeds

### ΔΔG Thresholds
| Threshold | Interpretation |
|-----------|----------------|
| ΔΔG ≥ 5 kcal/mol | High destabilization (structural collapse) |
| 2 ≤ ΔΔG < 5 kcal/mol | Moderate destabilization |
| ΔΔG < 2 kcal/mol | Minimal destabilization (functional mechanism) |

---

## Reproducibility

### Computational Environment

```
Platform: macOS 14.0 / Ubuntu 22.04
Python: 3.10+
FoldX: 5.1
Rosetta: 2023.49
```

### Runtime Estimates

| Analysis | Configuration | Time |
|----------|---------------|------|
| FoldX (96 variants × 3 seeds) | 10 replicates | ~15 hours |
| Rosetta (96 variants × 3 seeds) | 17 iterations, 10 cores | ~14 hours |

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

**Note**: FoldX and Rosetta require separate academic licenses.

---

## Contact

**Corresponding Author:**
Youngkyu Shim, M.D.
Department of Pediatrics, Korea University Ansan Hospital,
Korea University College of Medicine

**Issues**: Please open a GitHub issue for questions or bug reports

---

## Acknowledgments

- AlphaFold team at DeepMind for AlphaFold3
- FoldX team at CRG Barcelona
- RosettaCommons for Rosetta software
- ClinVar database for variant annotations
