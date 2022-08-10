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
    
    i = 0
    while reader:
        i+=1
        
        perc_comp = (1.0*i*chunk_size)/(L.mean())
        
        df_depth = df_depth_reader.get_chunk(chunk_size)
        df_refreq = df_refreq_reader.get_chunk(chunk_size)

        if df_depth.shape[0] < chunk_size:
            reader=False
            sys.stderr.write("Complete")
        else:
            sys.stderr.write("%s percent complete \n" % perc_comp)
        
        df_refreq = df_refreq.mask(df_depth < min_depth)   

        df_refcount = df_depth*df_refreq
        df_refreq_2 = df_refcount/(df_depth - 1)
       # df_refreq = df_refreq.where(df_refreq_2 > .05,0)
        #df_refreq = df_refreq.where(df_refreq_2 < .95,1)       
       # df_refreq_2 = df_refreq_2.where(df_refreq_2 > .05,0)
      #  df_refreq_2 = df_refreq_2.where(df_refreq_2 < .95,1)
        df_refreq_2 = df_refreq_2.mask(df_depth < min_depth)   

        df_altcount = df_depth*(1-df_refreq)
        df_altfreq = 1 - df_refreq
        df_altfreq_2 = df_altcount/(df_depth - 1)
     #   df_altfreq = df_altfreq.where(df_altfreq > .05,0)
     #   df_altfreq = df_altfreq.where(df_altfreq < .95,1)
     #   df_altfreq_2 = df_altfreq_2.where(df_altfreq_2 > .05,0)
     #   df_altfreq_2 = df_altfreq_2.where(df_altfreq_2 < .95,1)
        df_altfreq_2 = df_altfreq_2.mask(df_depth < min_depth)   

        for S in samples_tuples:
            
            s1 = S[0]
            s2 = S[1]
            s1s2 = (df_refreq[s1]*df_altfreq_2[s2] + df_refreq_2[s1]*df_altfreq[s2])

            pairwise_dic[(s1,s2)] += s1s2.sum()
            num_comps_dic[(s1,s2)] += s1s2.notna().sum()
    
    pi = {s:pairwise_dic[s]/num_comps_dic[s] for s in samples_tuples}
    
    pi_r = {(s[1],s[0]):pi[s] for s in samples_tuples}
    
    for s in pi_r.keys():
        pi[s] = pi_r[s]
    
    df = pd.Series(pi).unstack()
    
    dates = config.dates
    species_dates = dates[[d for d in dates.index if d in df.index]]
    
    df = df.loc[species_dates.index,species_dates.index]
    
    return df

if __name__ == "__main__":
    import sys
    
    species = sys.argv[1]
    pi = calculate_pi(species)
    if host != None:
        pi.to_csv("%s/pi/%s/%s/%s_pi.txt" % (analysis_dir,config.cohort,host,species))
    else:
        pi.to_csv("%s/pi/%s/%s_pi.txt" % (analysis_dir,config.cohort,species))
        