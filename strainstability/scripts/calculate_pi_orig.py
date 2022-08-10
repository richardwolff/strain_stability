import numpy as np
import pandas as pd
import config
import sys
import parse_midas_data
import itertools

host = config.host
analysis_dir = config.analysis_directory
min_depth = config.min_depth_pi_site
frac_sites = config.min_prev_pi_site
postproc_idx_dic = config.postproc_idx_dic
good_species = parse_midas_data.parse_good_species_list()

## implement pi estimate from Schloissnig (2013)
def calculate_pi(species): 
    
    snps_dir = "%ssnps/%s" % (config.data_directory,species)
    snps_summary = pd.read_csv("%s/snps_summary.txt" % snps_dir,sep="\t",index_col=0)
    L = snps_summary["covered_bases"]
   
    chunk_size = 40000
    samples_host = list(pd.read_csv("%s/snps_depth.txt.bz2" % snps_dir,sep="\t",index_col=0, nrows=0))
    samples_tuples = list(itertools.combinations(samples_host, 2)) + [(s1,s1) for s1 in samples_host]
    
    pi_vals = pd.DataFrame(0,index=samples_host,columns=samples_host)
    
    df_depth_reader = pd.read_csv("%s/snps_depth.txt.bz2" % snps_dir,sep="\t",index_col=0, iterator=True,low_memory=False)
    if config.MIDAS_ver == "1.3.2":
        df_refreq_reader = pd.read_csv("%s/snps_freq.txt.bz2" % snps_dir,sep="\t",index_col=0, iterator=True,low_memory=False)
     
    else:
        df_refreq_reader = pd.read_csv("%s/snps_ref_freq.txt.bz2" % snps_dir,sep="\t",index_col=0, iterator=True,low_memory=False)

    df_depth = df_depth_reader.get_chunk(0)
    df_refreq = df_refreq_reader.get_chunk(0)

    reader=True
   
    num_samples = len(samples_host)
    pairwise_dic = {s:0 for s in samples_tuples}
    num_comps_dic = {s:0 for s in samples_tuples}
    
    var_sites_df = pd.DataFrame(columns=samples_host)
    var_sites_depth_df = pd.DataFrame(columns=samples_host)
    i = 0
    
    G=pd.Series(0,index=samples_host)

    while reader:

        df_depth = df_depth_reader.get_chunk(chunk_size)
        df_refreq = df_refreq_reader.get_chunk(chunk_size)

        if df_depth.shape[0] < chunk_size:
            reader=False
            print("Complete")

        medians = df_depth.T.median()
        good_inds = medians[medians >= min_depth].index
        df_depth = df_depth.loc[good_inds]
        df_refreq = df_refreq.loc[good_inds]

        df_refreq = df_refreq*((1*(df_depth > min_depth)).replace(0,np.nan))
        df_depth = df_depth*(1*(df_depth > min_depth)).replace(0,np.nan)
        G+=1*(df_depth > min_depth).sum()

        df_refreq = df_refreq*((1*(df_refreq >= .05)))
        df_refreq = df_refreq*((1*(df_refreq <= .95)))

        df_refcount = df_depth*df_refreq
        df_altcount = df_depth*(1-df_refreq)

        df_pimat = 2*df_refreq*(df_altcount/(df_depth - 1))
        pi_vals += df_pimat.sum()   

        #host_refreq = df_refreq[samples_host]
        #var_sites = np.logical_and(host_refreq > 0,host_refreq < 1).T.sum() > L*frac_sites
        #var_sites_df = var_sites_df.append(host_refreq.loc[var_sites])
        #var_sites_depth_df = var_sites_depth_df.append(df_depth[samples_host].loc[var_sites])

        #S += sum(1*(var_sites))
        i+=1
        print(pi_vals)
        sys.stderr.write("step %s complete \n" % i)
    #pi = (pi_vals/G).loc[tp_matched[samples_host].sort_values().index]
    pi = pi_vals/G
    return(pi)

if __name__ == "__main__":
    import sys
    
    species = sys.argv[1]
    pi = calculate_pi(species)
    if host != None:
        pi.to_csv("%s/pi/%s/%s/%s_pi.txt" % (analysis_dir,config.cohort,host,species))
    else:
        pi.to_csv("%s/pi/%s/%s_pi.txt" % (analysis_dir,config.cohort,species))
        