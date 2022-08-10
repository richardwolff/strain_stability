host_flag=''
db_state_flag=''
cohort_flag=''

print_usage() {
  printf "Usage: \n specify host with -h \n specific cohort with -c \n specify midas db state with -d \n \n "
}

while getopts 'h:d:c:u' flag; do
  case "${flag}" in
    s) host_flag="${OPTARG}" ;;
    b) db_state_flag="${OPTARG}" ;;
    f) cohort_flag="${OPTARG}" ;;
    u) print_usage
       exit 1 ;;
  esac
done

export DESIRED_HOST=$host_flag
export DB_STATE=$db_state_flag
export DATA_COHORT=$cohort_flag