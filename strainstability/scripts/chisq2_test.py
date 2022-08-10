import state_utils
import pandas as pd
import numpy as np
import time
from scipy.stats import gamma
import config
import slm_utils
from scipy.stats import chisquare
from parse_midas_data import parse_good_species_list
import sys

## expected value of SLM at time t, given initial value of x0 at time 0 
def E_t(K,x0,t,tau):
    return K/(1+((K-x0)/x0)*np.exp(-t/tau))

## calculates the chi-squared test statistic for the goodness of fit test
def calc_chisq(obs_data,strain_dates,dates_diff,train_num):
    
    ## fit the SLM on the first train_num observations
    params = slm_utils.fit_SLM_params(obs_data,n=train_num)

    ## initialize SLM 
    S = slm_utils.slm(sigma=params["sigma"],K=params["K"],tau=1,delta_t=1.0/100)

    ## form quin_bins number of bins, expect to have at least 5 (= 4 + 1) observations in each bin, per recommendations re: chisquare test
    ## https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.chisquare.html
    quin_bins = int(len(obs_data)/4)
    #quin_bins = 13
    
    bin_list = []
    
    ## get bin of initial value w/r/t stationary distribution
    init_bin = int(np.digitize(obs_data.iloc[train_num],S.afd.ppf(np.linspace(0,100,quin_bins + 2)/100)[1:])) - 1
    bin_list.append(init_bin)
    
    ## now get bins of all subsequent timepoints 
    ## essentially, build an empirical probability distribution for abundance at time t + delta_t given abundance at t
    ## then, figure out which quintile bin of this empirical distribution the actual observation lies in
    for i in range(train_num,len(obs_data) - 1):
        
        T = dates_diff[i+1]
        num_iters = int(1.0*T/(1.0*S.delta_t))

        init_val = obs_data.iloc[i]
        num_reps = 100
        S.set_init_val(obs_data.iloc[i])
        S.run_sim(num_iters,num_reps,record_steps=False)

        bin_list.append(np.digitize(obs_data.iloc[i+1],np.percentile(S.x_i,np.linspace(0,100,quin_bins))))
    
    ## calculate chisquare test statistic relative to uniform expectation
    return chisquare(np.unique(bin_list,return_counts=True)[1])[1]  

#if __name__ == "main":
opt = sys.argv[1]
data_dir = config.data_directory
dates = config.dates
analysis_dir = config.analysis_directory
host = config.host
spec_df = pd.read_csv("%s/species/relative_abundance.txt.bz2" % data_dir,index_col=0,sep="\t")
strain_df = pd.read_csv("strains_%s.csv" % host,index_col=0)

good_species = parse_good_species_list()

chisq_dic = {}

dates_diff = dates.diff()

if opt == "strain":
    for strain in strain_df.index[::-1]:
      #  try:
        species = strain[:-2]
        good_samples = pd.read_csv("%s/snps/%s/snps_summary.txt" % (data_dir,species),index_col=0,sep="\t").index
        sys.stderr.write("Processing %s \n" % species)
        obs_data = strain_df.loc[strain]
        #print(obs_data)
        #obs_data = obs_data[good_samples]
        obs_data = obs_data[obs_data.notna()]

        #if len(obs_data) > 20:
        strain_dates = dates.loc[obs_data.index].sort_values()
        
        
        ## be sure observed data is in correct chronological order                
        obs_data = obs_data.loc[strain_dates.index]
        strain_dates_diff = strain_dates.diff()
        #print(strain_dates_diff)
        train_num = int(.333*len(obs_data)) 
        #train_num = np.argmin(abs((strain_dates.values - strain_dates[-1]/3)))
        #chisq_dic[strain] = calc_chisq(obs_data,strain_dates_diff,train_num)

        chisq_dic[strain] = np.median([calc_chisq(obs_data,strain_dates,strain_dates_diff,train_num) for _ in range(1)])
        print(chisq_dic[strain])
        print((pd.Series(chisq_dic) > .05).mean())
        print("\n")
       # except:
       #     pass
        
if opt == "species":
    for species in good_species:
       # try:
        good_samples = pd.read_csv("%s/snps/%s/snps_summary.txt" % (data_dir,species),index_col=0,sep="\t").index
        sys.stderr.write("Processing %s \n" % species)
        obs_data = spec_df.loc[species]
        obs_data = obs_data[good_samples]
        obs_data = obs_data[obs_data.notna()]

      #  if len(obs_data) > 20:
        spec_dates = dates.loc[obs_data.index].sort_values()

        ## be sure observed data is in correct chronological order
        obs_data = obs_data.loc[spec_dates.index]
        spec_dates_diff = spec_dates.diff()
        #train_num = np.argmin(abs((spec_dates.values - spec_dates[-1]/2)))

        train_num = int(.333*len(obs_data)) 
        #chisq_dic[species] = calc_chisq(obs_data,spec_dates_diff,train_num)
        chisq_dic[species] = np.median([calc_chisq(obs_data,spec_dates,spec_dates_diff,train_num) for _ in range(1)])
        print(chisq_dic[species])
        print((pd.Series(chisq_dic) > .05).mean())
        print("\n")
       # except:
       #     pass
        
if opt == "strain":
    pd.Series(chisq_dic).to_csv("%s/chisq/%s/%s_strain_chisq_test_cross.txt" % (config.analysis_directory,config.cohort,config.host))
    
elif opt == "species":
    pd.Series(chisq_dic).to_csv("%s/chisq/%s/%s_species_chisq_test.txt" % (config.analysis_directory,config.cohort,config.host))