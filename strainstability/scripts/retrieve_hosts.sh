readarray host_arr<${HOME}/strain_stability_revisions/strainstability/metadata/host_lists/${DATA_COHORT}_hosts.txt
host_arr=( null ${host_arr[@]} )
echo ${host_arr[@]}
