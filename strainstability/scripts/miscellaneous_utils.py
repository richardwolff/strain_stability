import pandas as pd
import config  
    
def return_dates(host):    
    all_dates = pd.read_csv("%s/host_dates/Poyet_dates.txt" % config.metadata_directory,index_col=0)
    dates = all_dates[all_dates["subject"] == host].sort_values("timepoint")["timepoint"]   
    return dates

def return_host_samples(host):
    
    host_samples_file = "%s/host_samples/poyet_%s.txt" % (config.metadata_directory,host)
    f = open(host_samples_file, "r")
    host_samples = f.readlines()
    host_samples = [fl.strip() for fl in host_samples] 
    return host_samples