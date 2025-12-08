# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2024-12-08

### Added
- Initial release accompanying manuscript submission
- AlphaFold3 structure prediction pipeline
- FoldX 5.1 stability analysis workflow
- Rosetta cartesian_ddG analysis workflow
- Statistical analysis scripts (Mann-Whitney U, ROC, effect sizes)
- 3 AlphaFold3 structural seeds (PDB format)
- Supplementary Methods documentation

### Scripts
- `01_structure_analysis/` - AlphaFold3 output analysis and seed selection
- `02_foldx/` - FoldX pipeline and result parsing
- `03_rosetta/` - Rosetta cartesian_ddG pipeline and parsing
- `04_analysis/` - Result combination and statistical analysis
