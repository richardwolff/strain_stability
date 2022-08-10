## Utility for retrieving and passing necessary environmental state variables 
## 
## The set_*() functions in this module do not actually write environmental state variables, they only set 
## these variables for child processes.
##
## To set the state variables, use the set_env_state.sh script. 

import os
import sys

def get_host():
    
    if "DESIRED_HOST" in os.environ.keys():
        host = os.environ["DESIRED_HOST"]
    else:
        host = None
    
    return host

def set_host(host):
    
    os.environ["DESIRED_HOST"] = host
    
    return host

def remove_host(host):
    
    del os.environ["DESIRED_HOST"] 

def get_midas_db_type():
    
    if "DB_TYPE" in os.environ.keys():
        db_type = os.environ["DB_TYPE"]
    else:
        db_type = None
    
    return db_type    

def set_midas_db_type(db_type):
    
    os.environ["DB_TYPE"] = db_type
    
    return db_type

def remove_midas_db_type():
    
    del os.environ["DB_TYPE"] 

def get_cohort():
    
    if "DATA_COHORT" in os.environ.keys():
        cohort = os.environ["DATA_COHORT"]
    else:
        cohort = None
    
    return cohort    

def set_cohort(cohort):
    
    os.environ["DATA_COHORT"] = cohort
    
    return cohort

def remove_cohort():
    
    del os.environ["DATA_COHORT"] 

def get_state_vars():
    
    state_vars = {}
    state_vars["DATA_COHORT"] = get_cohort()
    state_vars["DB_TYPE"] = get_midas_db_type()
    state_vars["DESIRED_HOST"] = get_host()
    
    return state_vars

def get_MIDAS_DB_dir():
    import config
    return config.midas_directory

def get_output_dir():
    import config
    return config.data_directory

## retrieves a dictionary of global attributes from config, gives them names to be exported as bash variables
## e.g. get_config_atts("MIDAS_DB" = "midas_directory", "OUTPUT" = "data_directory")

def get_config_named_atts(**kwargs):
    import config
    
    var_dic = {}
    #for key, value in kwargs.items():
    #    print("{0},{1}".format(key, value))

    for key, value in kwargs.items():
        var_dic[key] = getattr(config,value)
        
    return var_dic
    
if __name__ == "__main__":
    
    host = get_host()
    db_type = get_midas_db_type()
    cohort = get_cohort()

    if cohort == None:
        raise CohortError("DATA_COHORT environmental variable not set")
    else:
        sys.stderr.write("Data cohort: %s \n" % cohort)    
    
    if host == None:
        raise HostError("DESIRED_HOST environmental variable not set")
    else:
        sys.stderr.write("Desired host: %s \n" % host)

    if db_type == None:
        raise dbTypeError("DB_TYPE environmental variable not set")
    else:
        sys.stderr.write("DB type: %s \n" % db_type)
   
