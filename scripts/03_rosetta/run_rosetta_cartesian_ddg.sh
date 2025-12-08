#!/bin/bash

################################################################################
# Rosetta Cartesian DDG Protocol - HECW2 Variant Analysis
################################################################################
# Evaluate point mutations using cartesian_ddg with local backbone flexibility
# 10 background structures × 5 iterations = 50 evaluations per mutation
################################################################################

set -e

BASE_DIR="/Users/ykshim2025/Desktop/Code2025/HECW2"
ROSETTA_DIR="/Users/ykshim2025/Desktop/Code2025/FBXO11_SKP1/foldx_working/rosetta"
ROSETTA_BIN="${ROSETTA_DIR}/source/bin/cartesian_ddg.macosclangrelease"
ROSETTA_DB="${ROSETTA_DIR}/database"

# Input/Output directories
BACKGROUND_DIR="${BASE_DIR}/analysis/rosetta/relaxed_backgrounds"
MUTATION_FILE="${BASE_DIR}/rosetta_mutations.txt"
OUTPUT_DIR="${BASE_DIR}/analysis/rosetta/cartesian_ddg"
mkdir -p "${OUTPUT_DIR}"

# Parameters
NUM_ITERATIONS=5
NUM_PARALLEL=10

echo "================================================================================"
echo "ROSETTA CARTESIAN_DDG - HECW2 VARIANT ANALYSIS"
echo "================================================================================"
echo ""
echo "Base directory: ${BASE_DIR}"
echo "Background structures: ${BACKGROUND_DIR}"
echo "Mutation file: ${MUTATION_FILE}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "Parameters:"
echo "  - Iterations per background: ${NUM_ITERATIONS}"
echo "  - Parallel jobs: ${NUM_PARALLEL}"
echo ""
echo "Protocol: cartesian_ddg with 6-residue backbone neighborhood and ref2015_cart"
echo "Runtime estimate: 6-12 hours with parallel processing"
echo ""

# Check Rosetta binary
if [ ! -f "${ROSETTA_BIN}" ]; then
    echo "❌ ERROR: Rosetta cartesian_ddg binary not found at ${ROSETTA_BIN}"
    exit 1
fi

echo "✅ Using cartesian_ddg: ${ROSETTA_BIN}"
echo ""

# Check mutation file
if [ ! -f "${MUTATION_FILE}" ]; then
    echo "❌ ERROR: Mutation file not found at ${MUTATION_FILE}"
    exit 1
fi

# Count mutations (exclude header and comments)
NUM_MUTATIONS=$(grep -v "^#" "${MUTATION_FILE}" | grep -v "^total" | grep -v "^$" | wc -l | tr -d ' ')
echo "✅ Found ${NUM_MUTATIONS} mutations in ${MUTATION_FILE}"
echo ""

# Check background structures
NUM_BACKGROUNDS=$(ls "${BACKGROUND_DIR}"/relaxed_bg*.pdb 2>/dev/null | wc -l | tr -d ' ')
if [ "${NUM_BACKGROUNDS}" -eq 0 ]; then
    echo "❌ ERROR: No relaxed background structures found in ${BACKGROUND_DIR}"
    echo "Please run run_rosetta_fastrelax.sh first"
    exit 1
fi

echo "✅ Found ${NUM_BACKGROUNDS} relaxed background structures"
echo ""

TOTAL_ANALYSES=$((NUM_MUTATIONS * NUM_BACKGROUNDS * NUM_ITERATIONS))
echo "Total analyses: ${NUM_MUTATIONS} mutations × ${NUM_BACKGROUNDS} backgrounds × ${NUM_ITERATIONS} iterations = ${TOTAL_ANALYSES}"
echo ""

echo "================================================================================"
echo "STEP 1: Parse mutation file and create individual mutation files"
echo "================================================================================"
echo ""

MUT_DIR="${OUTPUT_DIR}/mutation_files"
mkdir -p "${MUT_DIR}"

# Parse mutation file and create individual .mut files
tail -n +3 "${MUTATION_FILE}" | grep -v "^#" | grep -v "^$" | while read line; do
    POS=$(echo "$line" | awk '{print $1}')
    CHAIN=$(echo "$line" | awk '{print $2}')
    WT=$(echo "$line" | awk '{print $3}')
    MUT=$(echo "$line" | awk '{print $4}')

    # Convert 1-letter to 3-letter amino acid codes
    case $WT in
        A) WT3="ALA";;  C) WT3="CYS";;  D) WT3="ASP";;  E) WT3="GLU";;
        F) WT3="PHE";;  G) WT3="GLY";;  H) WT3="HIS";;  I) WT3="ILE";;
        K) WT3="LYS";;  L) WT3="LEU";;  M) WT3="MET";;  N) WT3="ASN";;
        P) WT3="PRO";;  Q) WT3="GLN";;  R) WT3="ARG";;  S) WT3="SER";;
        T) WT3="THR";;  V) WT3="VAL";;  W) WT3="TRP";;  Y) WT3="TYR";;
    esac

    case $MUT in
        A) MUT3="ALA";;  C) MUT3="CYS";;  D) MUT3="ASP";;  E) MUT3="GLU";;
        F) MUT3="PHE";;  G) MUT3="GLY";;  H) MUT3="HIS";;  I) MUT3="ILE";;
        K) MUT3="LYS";;  L) MUT3="LEU";;  M) MUT3="MET";;  N) MUT3="ASN";;
        P) MUT3="PRO";;  Q) MUT3="GLN";;  R) MUT3="ARG";;  S) MUT3="SER";;
        T) MUT3="THR";;  V) MUT3="VAL";;  W) MUT3="TRP";;  Y) MUT3="TYR";;
    esac

    MUTFILE="${MUT_DIR}/p.${WT3}${POS}${MUT3}.mut"

    # Create mutation file in Rosetta format
    cat > "${MUTFILE}" << EOF
total 1
1
${WT} ${POS} ${CHAIN} PIKAA ${MUT}
EOF

    echo "  Created: p.${WT3}${POS}${MUT3}.mut"
done

NUM_MUTFILES=$(ls "${MUT_DIR}"/*.mut 2>/dev/null | wc -l | tr -d ' ')
echo ""
echo "✅ Created ${NUM_MUTFILES} mutation files"
echo ""

echo "================================================================================"
echo "STEP 2: Generate cartesian_ddg jobs for all mutations and backgrounds"
echo "================================================================================"
echo ""

JOB_LIST="${OUTPUT_DIR}/cartesian_ddg_jobs.txt"
rm -f "${JOB_LIST}"

# For each mutation
for MUTFILE in "${MUT_DIR}"/*.mut; do
    MUTNAME=$(basename "${MUTFILE}" .mut)

    # For each background structure
    for BG_PDB in "${BACKGROUND_DIR}"/relaxed_bg*.pdb; do
        BG_NUM=$(basename "${BG_PDB}" .pdb | sed 's/relaxed_bg//')

        WORK_DIR="${OUTPUT_DIR}/results/${MUTNAME}_bg${BG_NUM}"
        mkdir -p "${WORK_DIR}"

        # Copy input files
        cp "${BG_PDB}" "${WORK_DIR}/input.pdb"
        cp "${MUTFILE}" "${WORK_DIR}/mutations.mut"

        # Create cartesian_ddg command
        # Using 6-residue backbone neighborhood, ref2015_cart, 5 iterations
        CMD="cd ${WORK_DIR} && ${ROSETTA_BIN} \
            -database ${ROSETTA_DB} \
            -s input.pdb \
            -ddg:mut_file mutations.mut \
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
            -mute all \
            > ${MUTNAME}_bg${BG_NUM}.log 2>&1"

        echo "${CMD}" >> "${JOB_LIST}"
    done
done

TOTAL_JOBS=$(wc -l < "${JOB_LIST}")
echo "Total jobs to run: ${TOTAL_JOBS}"
echo "Parallel processes: ${NUM_PARALLEL}"
echo ""

echo "================================================================================"
echo "STEP 3: Run cartesian_ddg with GNU parallel"
echo "================================================================================"
echo ""

# Check for GNU parallel
if command -v parallel &> /dev/null; then
    echo "Using GNU parallel for job execution..."
    cat "${JOB_LIST}" | parallel -j ${NUM_PARALLEL} --progress
    echo ""
    echo "✅ cartesian_ddg completed for all ${TOTAL_JOBS} analyses"
else
    echo "⚠️  GNU parallel not found, running sequentially (this will take a long time)..."
    while IFS= read -r cmd; do
        eval "${cmd}"
    done < "${JOB_LIST}"
    echo ""
    echo "✅ cartesian_ddg completed for all ${TOTAL_JOBS} analyses"
fi

echo ""
echo "================================================================================"
echo "STEP 4: Collect and summarize results"
echo "================================================================================"
echo ""

# Count successful runs
RESULT_FILES=$(find "${OUTPUT_DIR}/results" -name "ddg_predictions.out" 2>/dev/null | wc -l | tr -d ' ')
echo "Found ${RESULT_FILES} ddg_predictions.out files"

if [ "${RESULT_FILES}" -gt 0 ]; then
    echo ""
    echo "Sample results (first 5):"
    find "${OUTPUT_DIR}/results" -name "ddg_predictions.out" 2>/dev/null | head -5 | while read file; do
        echo ""
        echo "$(dirname $file):"
        head -5 "$file"
    done
fi

echo ""
echo "================================================================================"
echo "SUMMARY"
echo "================================================================================"
echo ""
echo "Successfully completed ${RESULT_FILES}/${TOTAL_JOBS} cartesian_ddg analyses"
echo "Output location: ${OUTPUT_DIR}/results/"
echo ""
echo "Next step: Parse results with scripts/parse_rosetta_ddg.py"
echo ""
