# Code for running MSA portion of the AlphaFold3 pipeline on Virtual Machine
# Written by Miles Woodcock-Girard for Drew Lab at UIC

#!/bin/bash -l

# Define base AF3 macro-directory
BASE_DIR="/home/exouser"

# Define AF3 root directory
AF3_DIR="$BASE_DIR/alphafold3"

# Define input, output directories
INPUT_DIR="$BASE_DIR/input_candida"
OUTPUT_DIR="$BASE_DIR/output_candida"

# Define AF3 parameter directory
PARAM_DIR="$BASE_DIR/af3_params"

# Define AF3 database directory
AF3_DATA="$BASE_DIR/public_databases"



#cd $BASE_DIR

# Get number of files in input directory
#NUM_FILES=$(ls -1 | wc -l)

#FILE_START=0
#FILE_END=$((FILE_START + NNODES - 1))

# Make array of files
#file_array=($(ls ${INPUT_DIR} | sort | sed -n "$((FILE_START + 1)),$((FILE_END + 1))p"))
file_array=($(ls ${INPUT_DIR} | sort))

run_af3() {
	INPUT_FILE="$1"
	INPUT_STEM="${INPUT_FILE%.*}"
	OUTPUT_STEM="${INPUT_STEM,,}"
	OUTPUT_CURR="$OUTPUT_DIR/$OUTPUT_STEM"
	OUTPUT_JSON="$OUTPUT_CURR/${OUTPUT_STEM}_data.json"
	LOG_FILE="${OUTPUT_STEM}_log.txt"

	if [ -f "$OUTPUT_JSON" ]; then
		echo "MSA found for $INPUT_FILE. Skipping ..."
		return
	fi

	echo "Processing $INPUT_FILE"

	docker run \
		--volume $INPUT_DIR:/root/input_af3 \
		--volume $OUTPUT_DIR:/root/output_af3 \
		--volume $PARAM_DIR:/root/models \
		--volume $AF3_DATA:/root/databases \
		alphafold3 \
		python3 run_alphafold.py \
		--norun_inference \
		--json_path=/root/input_af3/${INPUT_FILE} \
		--model_dir=/root/models \
		--db_dir=/root/databases \
		--output_dir=/root/output_af3 \
		>> $LOG_FILE 2>&1
}

export -f run_af3
export INPUT_DIR OUTPUT_DIR PARAM_DIR AF3_DATA


parallel -j 8 run_af3 ::: "${file_array[@]}"
