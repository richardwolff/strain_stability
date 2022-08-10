. /u/local/Modules/default/init/modules.sh
module load anaconda2

conda activate strain_stability_env
source set_env_state.sh -c Poyet -i ae -d standard
exec bash
source declare_config_named_atts.sh
exec bash
echo $DESIRED_HOST