# Code for running MSA portion of the AlphaFold3 pipeline on HPC systems using PBS Pro
# Written by Miles Woodcock-Girard for Drew Lab at UIC

#!/bin/bash -l
#PBS -l select=10:system=polaris
#PBS -l place=scatter
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:eagle
#PBS -q prod
#PBS -N dc2_test
#PBS -A DirectContacts2
#PBS -M mwoodc2@uic.edu
#PBS -m bae

# Define base AF3 macro-directory
BASE_DIR="/lus/eagle/projects/DirectContacts2/miles/af3"

# Define local scratch directory to store MSA database
LOCAL_SCRATCH="/local/scratch/af3_db"
LOCAL_SCRATCH_COMPRESSED="/local/scratch/af3_db.tar.gz"

# Define AF3 root directory
AF3_DIR="$BASE_DIR/alphafold3"

# Define input, output directories
INPUT_DIR="$BASE_DIR/input_af3/json_input"
OUTPUT_DIR="$BASE_DIR/output_af3"

# Define AF3 parameter directory
PARAM_DIR="$BASE_DIR/af3_params"

# Define AF3 database directory
AF3_DATA="$BASE_DIR/af3_db"
AF3_DATA_COMPRESSED="$BASE_DIR/af3_db.tar.gz"



# Enable GPU-MPI
export MPICH_GPU_SUPPORT_ENABLED=1

cd $BASE_DIR

# Load necessary modules for running containers
module use /soft/modulefiles
module load spack-pe-base/0.8.1
module load apptainer

# Set proxy for internet access
export HTTP_PROXY=http://proxy.alcf.anl.gov:3128
export HTTPS_PROXY=http://proxy.alcf.anl.gov:3128
export http_proxy=http://proxy.alcf.anl.gov:3128
export https_proxy=http://proxy.alcf.anl.gov:3128

#ADDITIONAL_PATH=/opt/cray/pe/pals/1.2.12/lib
#module load cray-mpich-abi

# Get number of files in input directory
NUM_FILES=$(ls -1 | wc -l)
NNODES=`wc -l < $PBS_NODEFILE`

FILE_START=0
FILE_END=$((FILE_START + NNODES - 1))

# Make array of files, based on number of nodes
file_array=($(ls ${INPUT_DIR} | sort | sed -n "$((FILE_START + 1)),$((FILE_END + 1))p"))
#file_array=($(ls "$INPUT_DIR" | sort))
#num_files=${#file_array[@]}

# Split PBS node file into separate host files containing addressable node names
split --lines=1 --numeric-suffixes=1 --suffix-length=3 $PBS_NODEFILE local_hostfile.



# Loop over distinct nodes, start an AF3 instance on each
i=0
for lh in local_hostfile*; do

	INPUT_FILE="${file_array[$i]}"
        INPUT_STEM="${INPUT_FILE%.*}"

        OUTPUT_STEM="${INPUT_STEM,,}"
        OUTPUT_CURR="$OUTPUT_DIR/$OUTPUT_STEM"
        OUTPUT_CIF="$OUTPUT_CURR/${OUTPUT_STEM}_model.cif"

	echo $OUTPUT_CURR
	echo $OUTPUT_CIF

        if [ -f $OUTPUT_CIF ]; then
                echo "Model file found for $INPUT_FILE. Skipping"
		i=$((i + 1))
                continue
        fi

	LOG_FILE="node_${i}_log.txt"
	echo "$INPUT_FILE"

	time mpiexec -n 1 --ppn 1 --hostfile ${lh} \
		apptainer exec \
		--nv \
		--bind $INPUT_DIR:/root/input_af3 \
		--bind $OUTPUT_DIR:/root/output_af3 \
		--bind $PARAM_DIR:/root/models \
		--bind $AF3_DATA:/root/databases \
		$AF3_DIR/alphafold3.sif \
		python3 $AF3_DIR/run_alphafold.py \
		--norun_inference \
		--json_path=/root/input_af3/${INPUT_FILE} \
		--model_dir=/root/models \
		--db_dir=/root/databases \
		--output_dir=/root/output_af3 \
		>> "$LOG_FILE" 2>&1 &
	i=$((i + 1))

	sleep 10s
done

wait
