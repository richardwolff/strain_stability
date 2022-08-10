import state_utils

if __name__ == '__main__':
    
    ## TO-DO: fully modularize this so arbitrary inputs. Good enough for now. 
    output=state_utils.get_config_named_atts(MIDAS_DB="midas_directory",OUTDIR="data_directory",
                                             HOST_SAMPLES_FILE="host_samples_file",MIDAS_FILES_DIR="midas_files_directory",
                                            SF_DIR="strainfinder_directory")
    output_str_list = ["['%s']=%s" % (key,value) for key,value in output.items()]
    
    output_str_list.insert(0,"(")
    output_str_list.append(")")
    
    output_str = " ".join(output_str_list)
    print output_str
