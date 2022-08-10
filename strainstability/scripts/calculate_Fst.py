import numpy as np
import pandas as pd
import config
import time
import sys

cohort = config.cohort
host = config.host
data_dir = config.data_directory
analysis_dir = config.analysis_directory
species = sys.argv[1]
snps_dir = "%ssnps/%s" % (config.data_directory,species)

## implement pi within estimate from Schloissnig (2013)
chunk_size = 40000

samples_host = list(pd.read_csv("%s/snps_depth.txt.bz2" % snps_dir,sep="\t",index_col=0, nrows=0))

pi_vals = pd.DataFrame(0,index=samples_host,columns=samples_host)
df_depth_reader = pd.read_csv("%s/snps_depth.txt.bz2" % snps_dir,sep="\t",index_col=0, iterator=True,low_memory=False)
df_refreq_reader = pd.read_csv("%s/snps_ref_freq.txt.bz2" % snps_dir,sep="\t",index_col=0, iterator=True,low_memory=False)

df_depth = df_depth_reader.get_chunk(1)
df_refreq = df_refreq_reader.get_chunk(1)

t1 = time.time()
reader=True
i=0
G=pd.DataFrame(0,index=samples_host,columns=samples_host)

min_depth = 5

while reader:
    
    df_depth = df_depth_reader.get_chunk(chunk_size)
    df_refreq = df_refreq_reader.get_chunk(chunk_size)
    
    df_depth = df_depth[samples_host]
    df_refreq = df_refreq[samples_host]
    
    ## change the depth of all sites which have less than min_depth = 5
    ## to have an observed depth of 2
    ## The value of 2 is chosen so that the coverage correction will not have any 
    df_depth = df_depth.where(df_depth > min_depth,2)
    
    if df_depth.shape[0] < chunk_size:
        reader=False
        print("Complete")
    
    ## following Schloissnig, remove low freqeuncy variants
    df_refreq = df_refreq.where(df_refreq > .01,0)
    df_refreq = df_refreq.where(df_refreq < .99,1)
    ## 0 out sites with insufficient depth
    df_refreq = df_refreq.where(df_depth > min_depth, 0)   
    df_refcount = df_refreq*df_depth
    
    df_altfreq = 1 - df_refreq
    ## 0 out sites with insufficient depth
    df_altfreq = df_altfreq.where(df_depth > min_depth, 0)   
    df_altcount = df_altfreq*df_depth
    
    ## coverage correction
    ## as ref, alt, ref2, and alt2 all have freq 0, these sites will not contribute to the final calculation of pi
    ## as the denominator (G) only includes sites with depth > min_depth
    
    df_altfreq2 = df_altcount/(df_depth - 1)
    ## 0 out sites with insufficient depth
    df_altfreq2 = df_altfreq2.where(df_depth > min_depth, 0) 
    df_refreq2 = df_refcount/(df_depth - 1)
    ## 0 out sites with insufficient depth
    df_refreq2 = df_refreq2.where(df_depth > min_depth, 0) 
    
    df_pi_mat = (df_refreq.T @ df_altfreq2) + (df_altfreq.T @ df_refreq2)
    
    #df_pi_mat = (df_refreq.T @ df_altfreq) + (df_altfreq.T @ df_refreq)
    
    ## count number of sites where each comparison of samples has sufficient depth
    G += (1*(df_depth.T > min_depth)) @ (1*(df_depth > min_depth))
    
    pi_vals += df_pi_mat

    i+=1
    
    t2 = time.time()
    print(pi_vals.head())
    #sys.stderr.write(f"step {i}, time elapsed {t2 - t1}. Total sites: {G.iloc[0][0]/(chunk_size*i)} \n")
    
pi = pi_vals/G

print(pi.head())

pi.to_csv("%s/pi/%s/%s/%s_pi.txt" % (analysis_dir,cohort,host,species))

sys.stderr.write("Files written \n \n \n")
