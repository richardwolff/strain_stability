H=$1

module load singularity 

source ./set_env_state.sh -c Poyet -i $H -d standard

#hosts=( $( ./retrieve_hosts.sh ) )

#export DESIRED_HOST=${hosts[$H]}

source ./declare_config_named_atts.sh

mkdir -p $OUTDIR/snps
mkdir -p $OUTDIR/genes

readarray -t samples < $HOST_SAMPLES_FILE

input_host_files=( "${samples[@]/#/${MIDAS_FILES_DIR}}" )

input_host_list=$(IFS=,; echo "${input_host_files[*]}")

singularity exec $H2_CONTAINER_LOC/MIDAS-mod.sif merge_midas.py snps $OUTDIR/snps -i $input_host_list -t list --sample_depth 5 --site_depth 3 --min_samples 25 --max_species 150 --site_prev 0.0 --threads 10

#> $OUTDIR/snps/species.log

singularity exec $H2_CONTAINER_LOC/MIDAS-mod.sif merge_midas.py genes $OUTDIR/genes -i $input_host_list -t list --sample_depth 5 --min_samples 25

#  > $OUTDIR/genes/species.log





