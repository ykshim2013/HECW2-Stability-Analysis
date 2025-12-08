#!/bin/bash

#############################################################################
# FoldX Pipeline for HECW2 Variant Analysis
#
# Purpose: Automated FoldX BuildModel + Stability analysis for all 3 seeds
#
# Usage:
#   ./run_foldx_pipeline.sh test      # Run with 3 test mutations
#   ./run_foldx_pipeline.sh full      # Run all 31 variants
#   ./run_foldx_pipeline.sh seed1     # Run only seed1
#   ./run_foldx_pipeline.sh seed2     # Run only seed2
#   ./run_foldx_pipeline.sh seed3     # Run only seed3
#
# Timeline:
#   Test mode: ~5-10 minutes (3 mutations × 10 runs)
#   Full mode: ~9-12 hours (31 mutations × 10 runs × 3 seeds)
#
# Output:
#   analysis/foldx/seed{1,2,3}_foldx/
#     ├── WT_*.pdb (10 wildtype models)
#     ├── {variant}_*.pdb (310 mutant models per seed)
#     ├── Average_*.fxout (stability analysis results)
#     ├── Dif_*.fxout (ΔΔG values)
#     └── Raw_*.fxout (detailed energies)
#
#############################################################################

set -euo pipefail

# Configuration
FOLDX="/Applications/YASARA.app/Contents/yasara/plg/foldx"
PROJECT_DIR="/Users/ykshim2025/Desktop/Code2025/HECW2"
FOLDX_DIR="${PROJECT_DIR}/analysis/foldx"
N_RUNS=10

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Verify FoldX executable
if [ ! -f "$FOLDX" ]; then
    log_error "FoldX not found at: $FOLDX"
    exit 1
fi

log_info "Using FoldX: $FOLDX"
$FOLDX --version 2>&1 | head -3

# Parse command line arguments
MODE="${1:-test}"

case "$MODE" in
    test)
        SEEDS=("seed1")
        MUTATION_LIST="test_mutations.txt"
        log_info "TEST MODE: Running 3 mutations on seed1 only"
        ;;
    full)
        SEEDS=("seed1" "seed2" "seed3")
        MUTATION_LIST="individual_list_corrected.txt"
        log_info "FULL MODE: Running all 31 variants on all 3 seeds"
        ;;
    seed1|seed2|seed3)
        SEEDS=("$MODE")
        MUTATION_LIST="individual_list_corrected.txt"
        log_info "SINGLE SEED MODE: Running all 31 variants on $MODE"
        ;;
    *)
        log_error "Invalid mode: $MODE"
        echo "Usage: $0 {test|full|seed1|seed2|seed3}"
        exit 1
        ;;
esac

# Verify mutation list exists
if [ ! -f "${FOLDX_DIR}/${MUTATION_LIST}" ]; then
    log_error "Mutation list not found: ${FOLDX_DIR}/${MUTATION_LIST}"
    exit 1
fi

N_VARIANTS=$(wc -l < "${FOLDX_DIR}/${MUTATION_LIST}")
log_info "Mutation list: ${MUTATION_LIST} (${N_VARIANTS} variants)"

#############################################################################
# Function: Run FoldX BuildModel
#############################################################################
run_buildmodel() {
    local seed=$1
    local pdb_file=$2
    local output_dir=$3

    log_info "Starting BuildModel for ${seed}..."
    log_info "  Input PDB: ${pdb_file}"
    log_info "  Output dir: ${output_dir}"

    # Create output directory
    mkdir -p "$output_dir"

    # Copy required files to output directory
    cp "$pdb_file" "$output_dir/"
    cp "${FOLDX_DIR}/${MUTATION_LIST}" "$output_dir/individual_list.txt"

    # Get PDB filename
    pdb_name=$(basename "$pdb_file")

    # Change to output directory (FoldX outputs to current directory)
    cd "$output_dir"

    # Run BuildModel
    log_info "Running BuildModel (${N_RUNS} runs)..."
    $FOLDX --command=BuildModel \
           --pdb="$pdb_name" \
           --mutant-file=individual_list.txt \
           --numberOfRuns=$N_RUNS \
           --output-dir=./ \
           > buildmodel.log 2>&1

    # Check if models were generated
    n_models=$(ls -1 *.pdb 2>/dev/null | grep -v "$pdb_name" | wc -l)
    expected_models=$((N_VARIANTS * N_RUNS))

    if [ "$n_models" -eq "$expected_models" ]; then
        log_success "BuildModel complete: ${n_models}/${expected_models} models generated"
    else
        log_warning "BuildModel incomplete: ${n_models}/${expected_models} models generated"
        log_warning "Check buildmodel.log for details"
    fi

    cd "$PROJECT_DIR"
}

#############################################################################
# Function: Run FoldX Stability
#############################################################################
run_stability() {
    local seed=$1
    local output_dir=$2

    log_info "Starting Stability analysis for ${seed}..."

    cd "$output_dir"

    # Get list of mutant PDB files (exclude wildtype and original repaired)
    mutant_pdbs=$(ls -1 *.pdb 2>/dev/null | grep -v "Repair.pdb" | grep -v "^WT_")
    n_mutants=$(echo "$mutant_pdbs" | wc -l)

    log_info "Found ${n_mutants} mutant structures"

    # Run Stability on each mutant
    counter=0
    for pdb in $mutant_pdbs; do
        counter=$((counter + 1))

        if [ $((counter % 50)) -eq 0 ]; then
            log_info "Progress: ${counter}/${n_mutants} mutants analyzed"
        fi

        $FOLDX --command=Stability \
               --pdb="$pdb" \
               --output-dir=./ \
               >> stability.log 2>&1
    done

    log_success "Stability analysis complete: ${n_mutants} mutants analyzed"

    # Check for output files
    n_fxout=$(ls -1 *.fxout 2>/dev/null | wc -l)
    log_info "Generated ${n_fxout} .fxout files"

    cd "$PROJECT_DIR"
}

#############################################################################
# Main Pipeline
#############################################################################

log_info "========================================="
log_info "FoldX Pipeline Starting"
log_info "========================================="
log_info "Mode: $MODE"
log_info "Seeds: ${SEEDS[*]}"
log_info "Variants: $N_VARIANTS"
log_info "Runs per variant: $N_RUNS"
log_info "Expected models per seed: $((N_VARIANTS * N_RUNS))"
log_info "========================================="

# Record start time
start_time=$(date +%s)

# Process each seed
for seed in "${SEEDS[@]}"; do
    log_info ""
    log_info "========================================="
    log_info "Processing ${seed}"
    log_info "========================================="

    # Define paths
    repaired_pdb="${FOLDX_DIR}/repaired_wt/wt_${seed}.pdb"
    output_dir="${FOLDX_DIR}/${seed}_foldx"

    # Check if repaired PDB exists
    if [ ! -f "$repaired_pdb" ]; then
        log_error "Repaired PDB not found: $repaired_pdb"
        log_error "Please run RepairPDB first"
        exit 1
    fi

    # Step 1: BuildModel
    run_buildmodel "$seed" "$repaired_pdb" "$output_dir"

    # Step 2: Stability
    run_stability "$seed" "$output_dir"

    log_success "${seed} complete!"
done

# Calculate elapsed time
end_time=$(date +%s)
elapsed=$((end_time - start_time))
hours=$((elapsed / 3600))
minutes=$(((elapsed % 3600) / 60))
seconds=$((elapsed % 60))

log_info ""
log_info "========================================="
log_success "FoldX Pipeline Complete!"
log_info "========================================="
log_info "Total time: ${hours}h ${minutes}m ${seconds}s"
log_info ""
log_info "Next steps:"
log_info "  1. Check logs in analysis/foldx/{seed}_foldx/*.log"
log_info "  2. Parse results: python3 scripts/parse_foldx_results.py"
log_info "  3. Combine with Rosetta data for dual-method analysis"
log_info "========================================="
