cohort_flag=''
host_flag=''
db_type_flag=''

Help() {
   echo 
   echo "Set the DESIRED_HOST, DB_TYPE, and DATA_COHORT environmental variables here"
   echo
   echo "Syntax: . set_env_state.sh [-c|i|d|h]"
   echo "options:"
   echo "c     Cohort to process: Poyet, Korpela, Suez, Roodgar, or HMP. Env var: DATA_COHORT."
   echo "i     Individual host to process (e.g. am (Poyet), alien (Korpela), 10 (Suez)). Env var: DESIRED_HOST"
   echo "d     MIDAS database alignment type: standard or isolate. Env var: DB_TYPE"
   echo "h     Print this Help."
   echo
}

while getopts "c:i:d:h" flag; do 
    case $flag in   
        c) cohort_flag="${OPTARG}";
           echo "Data cohort: " $cohort_flag;
           export DATA_COHORT=$cohort_flag;;  
        i) host_flag="${OPTARG}";
           echo "Desired host: " $host_flag;
           export DESIRED_HOST=$host_flag;; 
        d) db_type_flag="${OPTARG}";
           echo "DB type: " $db_type_flag;
           export DB_TYPE=$db_type_flag;;
        h) Help;;
    esac
done


