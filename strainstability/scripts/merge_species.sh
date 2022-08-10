H=$1

module load singularity 

source ./set_env_state.sh -c Poyet -i $H -d standard

source ./declare_config_named_atts.sh

mkdir -p $OUTDIR/species

readarray -t samples < $HOST_SAMPLES_FILE

input_host_files=( "${samples[@]/#/${MIDAS_FILES_DIR}}" )

input_host_list=$(IFS=,; echo "${input_host_files[*]}")

singularity exec $H2_CONTAINER_LOC/MIDAS-mod.sif merge_midas.py species $OUTDIR/species -i $input_host_list -t list > $OUTDIR/species/species.log
