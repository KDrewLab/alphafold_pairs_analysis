# Code for running inference stage of the AlphaFold3 pipeline on Virtual Machine
# Written by Miles Woodcock-Girard for Drew Lab at UIC

#!/bin/bash -l

# Define base AF3 macro-directory
BASE_DIR="/home/exouser"

# Define AF3 root directory
AF3_DIR="$BASE_DIR/alphafold3"

# Define input, output directories
INPUT_DIR="$BASE_DIR/input_af3"
OUTPUT_DIR="$BASE_DIR/output_af3"

# Define AF3 model weights directory
PARAM_DIR="$BASE_DIR/af3_params"

# Make array of input files
file_array=($(ls ${INPUT_DIR} | sort))

# Define function for calling AlphaFold3 docker image
run_af3() {
	INPUT_FILE="$1"
	INPUT_STEM="${INPUT_FILE%.*}"
	INPUT_JSON="$INPUT_FILE/${INPUT_STEM}_data.json"
	OUTPUT_STEM="${INPUT_STEM,,}"
	OUTPUT_CURR="$OUTPUT_DIR/$OUTPUT_STEM"
	OUTPUT_CIF="$OUTPUT_CURR/${OUTPUT_STEM}_model.cif"
	LOG_FILE="${OUTPUT_STEM}_log.txt"

	if [ ! -f "$INPUT_DIR/$INPUT_JSON" ]; then
		echo "No AlphaFold3 data JSON found in $INPUT_FILE. Skipping ..."
		return
	fi

	if [ -f "$OUTPUT_CIF" ]; then
		echo "Models found for $INPUT_FILE. Skipping ..."
		return
	fi

	echo "Processing $INPUT_JSON"

	docker run \
		--volume $INPUT_DIR:/root/input_af3 \
		--volume $OUTPUT_DIR:/root/output_af3 \
		--volume $PARAM_DIR:/root/models \
		--gpus all \
		alphafold3 \
		python3 run_alphafold.py \
		--norun_data_pipeline \
		--json_path=/root/input_af3/${INPUT_JSON} \
		--model_dir=/root/models \
		--output_dir=/root/output_af3 \
		>> $LOG_FILE 2>&1
}

export -f run_af3
export INPUT_DIR OUTPUT_DIR PARAM_DIR AF3_DATA


# Loop through all directories in input folder, running AlphaFold3 on each
for file in "${file_array[@]}"; do
	run_af3 "$file"
done
