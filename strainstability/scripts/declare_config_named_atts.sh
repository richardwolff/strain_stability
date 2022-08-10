named_atts=$( python print_config_named_atts.py )

declare -A varnames=$named_atts

for name in "${!varnames[@]}"; 
    do export $name=${varnames[$name]}
done
