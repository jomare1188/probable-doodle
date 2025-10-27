#!/bin/bash
# HMMER search for transcription factors in sugarcane proteins
# Optimized for parallel execution with GNU parallel

set -e  # Exit on error

# Configuration
PROTEIN_FILE="/home/diegoj/rnaseq_diatraea/reference_genomes/sugarcane/annotation/SofficinarumxspontaneumR570_771_v2.1.protein.fa.gz"
HMM_DB="../db/TF.db.hmm"  # or OwnPlnTFDB.hmm.gz | from https://doi.org/10.1186/1471-2105-8-42
OUTPUT_DIR="../results/sugarcane"
TEMP_DIR="../temp_chunks"
NUM_CHUNKS=100  # Split protein file into 100 chunks for parallel processing
PARALLEL_JOBS=100  # Run 100 jobs in parallel (adjust based on your preference)
THREADS_PER_JOB=2  # Threads per hmmsearch job (100 jobs × 2 threads = 200 cores)

echo "========================================"
echo "HMMER Parallel Search - Sugarcane TFs"
echo "========================================"
echo "Configuration:"
echo "  Chunks: ${NUM_CHUNKS}"
echo "  Parallel jobs: ${PARALLEL_JOBS}"
echo "  Threads per job: ${THREADS_PER_JOB}"
echo "  Total cores used: $((PARALLEL_JOBS * THREADS_PER_JOB))"
echo "========================================"

# Create directories
mkdir -p ${OUTPUT_DIR}/parts
mkdir -p ${TEMP_DIR}
mkdir -p code/logs

# Check files
if [ ! -f ${PROTEIN_FILE} ]; then
    echo "Error: Protein file not found: ${PROTEIN_FILE}"
    exit 1
fi

if [ ! -f ${HMM_DB} ] && [ ! -f ${HMM_DB}.gz ]; then
    echo "Error: HMM database not found: ${HMM_DB}"
    exit 1
fi

# Prepare HMM database
if [[ ${HMM_DB} == *.gz ]] && [ -f ${HMM_DB} ]; then
    echo "Uncompressing HMM database..."
    HMM_DB_UNCOMPRESSED="${HMM_DB%.gz}"
    if [ ! -f ${HMM_DB_UNCOMPRESSED} ]; then
        gunzip -c ${HMM_DB} > ${HMM_DB_UNCOMPRESSED}
    fi
    HMM_DB=${HMM_DB_UNCOMPRESSED}
elif [ -f ${HMM_DB}.gz ]; then
    echo "Uncompressing HMM database..."
    gunzip -c ${HMM_DB}.gz > ${HMM_DB}
fi

# Press HMM database if needed
if [ ! -f ${HMM_DB}.h3m ]; then
    echo "Pressing HMM database..."
    hmmpress ${HMM_DB}
fi

# Prepare protein file
echo "Preparing protein file..."
if [[ ${PROTEIN_FILE} == *.gz ]]; then
    PROTEIN_UNCOMPRESSED="${TEMP_DIR}/proteins.fa"
    if [ ! -f ${PROTEIN_UNCOMPRESSED} ]; then
        echo "Uncompressing protein file..."
        gunzip -c ${PROTEIN_FILE} > ${PROTEIN_UNCOMPRESSED}
    fi
    PROTEIN_FILE=${PROTEIN_UNCOMPRESSED}
fi

# Count sequences
echo "Counting sequences..."
NUM_SEQS=$(grep -c "^>" ${PROTEIN_FILE})
echo "Total sequences: ${NUM_SEQS}"

# Split protein file into chunks
echo "Splitting protein file into ${NUM_CHUNKS} chunks..."
if [ ! -f ${TEMP_DIR}/chunk_001.fa ]; then
    seqkit split -p ${NUM_CHUNKS} -O ${TEMP_DIR} ${PROTEIN_FILE}
    # Rename files for easier processing
    cd ${TEMP_DIR}
    for f in proteins.part_*.fa; do
        num=$(echo $f | sed 's/proteins.part_//' | sed 's/.fa//')
        mv $f chunk_$(printf "%03d" $num).fa
    done
    cd - > /dev/null
fi

# Count chunks created
ACTUAL_CHUNKS=$(ls ${TEMP_DIR}/chunk_*.fa 2>/dev/null | wc -l)
echo "Created ${ACTUAL_CHUNKS} chunks"

if [ ${ACTUAL_CHUNKS} -eq 0 ]; then
    echo "Error: No chunks were created!"
    exit 1
fi

# Create GNU parallel command
echo "Running HMMER searches in parallel..."
echo "Start time: $(date)"

# Get absolute paths to avoid issues with parallel
ABS_HMM_DB=$(readlink -f ${HMM_DB})
ABS_OUTPUT_DIR=$(readlink -f ${OUTPUT_DIR})
ABS_TEMP_DIR=$(readlink -f ${TEMP_DIR})

# Run parallel searches with inline command
ls ${ABS_TEMP_DIR}/chunk_*.fa | \
    parallel -j ${PARALLEL_JOBS} --bar --eta \
    'chunk="{}"; \
    chunk_name=$(basename "$chunk" .fa); \
    echo "Processing $chunk_name..."; \
    hmmsearch \
        --cpu '"${THREADS_PER_JOB}"' \
        --domtblout "'"${ABS_OUTPUT_DIR}"'/parts/domtbl.$chunk_name.out" \
        --cut_ga \
        -Z '"${NUM_SEQS}"' \
        -o "'"${ABS_OUTPUT_DIR}"'/parts/all.$chunk_name.out" \
        "'"${ABS_HMM_DB}"'" \
        "$chunk" > "'"${ABS_OUTPUT_DIR}"'/parts/$chunk_name.log" 2>&1 || \
        echo "ERROR: Failed processing $chunk_name" >> "'"${ABS_OUTPUT_DIR}"'/parts/errors.log"'

echo "Parallel searches completed: $(date)"

# Check for errors
if [ -f ${OUTPUT_DIR}/parts/errors.log ]; then
    echo "WARNING: Some jobs failed. Check ${OUTPUT_DIR}/parts/errors.log"
    cat ${OUTPUT_DIR}/parts/errors.log
fi

# Verify output files exist
OUTPUT_COUNT=$(ls ${OUTPUT_DIR}/parts/domtbl.chunk_*.out 2>/dev/null | wc -l)
echo "Found ${OUTPUT_COUNT} output files from ${ACTUAL_CHUNKS} chunks"

if [ ${OUTPUT_COUNT} -eq 0 ]; then
    echo "ERROR: No output files were created!"
    echo "Check log files in ${OUTPUT_DIR}/parts/"
    exit 1
fi

# Concatenate results
echo "Concatenating results..."

# Concatenate domain tables (skip headers except first)
echo "Merging domain tables..."
first=1
for f in ${OUTPUT_DIR}/parts/domtbl.chunk_*.out; do
    if [ ! -f "$f" ]; then
        echo "Warning: Expected file not found: $f"
        continue
    fi
    
    if [ $first -eq 1 ]; then
        cat "$f" > ${OUTPUT_DIR}/domtbl_sugarcane.out
        first=0
    else
        grep -v "^#" "$f" >> ${OUTPUT_DIR}/domtbl_sugarcane.out
    fi
done

# Concatenate full outputs
echo "Merging full outputs..."
cat ${OUTPUT_DIR}/parts/all.chunk_*.out > ${OUTPUT_DIR}/all_sugarcane.out 2>/dev/null || \
    echo "Warning: Could not merge all outputs"

# Generate summary
echo "========================================"
echo "Analysis Complete!"
echo "========================================"

if [ -f ${OUTPUT_DIR}/domtbl_sugarcane.out ]; then
    HITS=$(grep -v "^#" ${OUTPUT_DIR}/domtbl_sugarcane.out | wc -l)
    UNIQUE_PROTEINS=$(grep -v "^#" ${OUTPUT_DIR}/domtbl_sugarcane.out | awk '{print $1}' | sort -u | wc -l)

    echo "Results:"
    echo "  Domain hits: ${HITS}"
    echo "  Unique proteins with TF domains: ${UNIQUE_PROTEINS}"
    echo "  Output files:"
    echo "    - ${OUTPUT_DIR}/domtbl_sugarcane.out"
    echo "    - ${OUTPUT_DIR}/all_sugarcane.out"
    echo ""

    # Extract TF IDs
    echo "Extracting TF identifiers..."
    grep -v "^#" ${OUTPUT_DIR}/domtbl_sugarcane.out | awk '{print $1}' | sort -u > ${OUTPUT_DIR}/TF_hits.ids
    echo "  - ${OUTPUT_DIR}/TF_hits.ids (${UNIQUE_PROTEINS} IDs)"
else
    echo "ERROR: Final output file was not created!"
    exit 1
fi

# Optional: Clean up temp files
read -p "Remove temporary chunk files? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Cleaning up temporary files..."
    rm -rf ${TEMP_DIR}
    rm -rf ${OUTPUT_DIR}/parts
    echo "Cleanup complete!"
fi

echo "========================================"
echo "Finished: $(date)"
echo "========================================"
