# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-XX-XX

### Added
- Initial release accompanying manuscript submission
- AlphaFold3 structure prediction pipeline
- FoldX 5.1 stability analysis workflow
- Rosetta cartesian_ddG analysis workflow
- Statistical analysis scripts (Mann-Whitney U, ROC, effect sizes)
- Figure generation scripts
- Complete dataset of 96 HECW2 variants
- Supplementary Methods documentation

### Data
- 96 HECW2 missense variants analyzed
  - 34 pathogenic/likely pathogenic
  - 22 benign/likely benign
  - 40 variants of uncertain significance (VUS)
- 3 AlphaFold3 structural seeds
- FoldX ΔΔG values (median ± IQR)
- Rosetta ΔΔG values (mean ± SD)

### Analysis
- Dual-method ΔΔG comparison
- Mechanism stratification (high/moderate/non-destabilization)
- Clinical phenotype correlation
- VUS reclassification candidates identified

---

## [Unreleased]

### Planned
- Integration with Missense3D structural analysis
- MD simulation validation for selected variants
- Web interface for variant lookup
