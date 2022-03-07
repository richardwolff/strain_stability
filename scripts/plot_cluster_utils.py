## python3
import pandas as pd
import numpy as np
import pickle
import config
import matplotlib.pyplot as plt
plt.style.use('ggplot')
from matplotlib import rc
rc('text', usetex=True)

SMALL_SIZE=15
MEDIUM_SIZE=20
rc('legend', fontsize=MEDIUM_SIZE)
rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
rc('ytick', labelsize=SMALL_SIZE) 
#colors = 10*['#348ABD', '#A60628', '#7A68A6', '#467821']
    
def get_clusters_snv_trajectories(snp_map):
    
    clusters = []
    for key in snp_map.keys():
        clusters.append(snp_map[key][1])
    clusters = list(set(clusters))
    
    freqs = {cluster:[] for cluster in clusters}
    for key in snp_map.keys():
        freqs[list(snp_map[key])[1]].append(np.array(list(snp_map[key])[2]/list(snp_map[key][3])))
    
    for key in freqs.keys():
        freqs[key] = np.array(freqs[key])
        
    return freqs

def plot_polarized_snv_trajectories(clusters_snv_trajectories,good_inds,ax,num_snps=1000,alpha=.01,):
    #fig,ax = plt.subplots(figsize=(12,8))
    
    for cluster in clusters_snv_trajectories.keys():
        cluster_snvs = clusters_snv_trajectories[cluster].T
        ax.plot(cluster_snvs[good_inds,:1000],alpha=alpha,color="grey");
    
    ax.set_xlabel("Timepoint")
    ax.set_ylabel("Frequency") 

def return_snp_strain_freqs(species,host):

    import pandas as pd
    
    output_directory = "/u/scratch/r/rwolff/strainfinder_input/Poyet/%s" % host
    filename_prefix = "%s/%s" % (output_directory, species)
    snp_locations = pickle.load(open(filename_prefix+".strainfinder.locations.p",'rb'))
    snp_alignment = pd.read_pickle(filename_prefix+".strainfinder.p")
    snp_samples = pickle.load(open(filename_prefix+".strainfinder.samples.p",'rb'))
    snp_samples = [elem.decode("utf-8") for elem in snp_samples]
    
    
    Poyet_samples = config.Poyet_samples 
    
    good_samples = Poyet_samples[host]
        
    outp=pd.read_pickle("%sclusters/%s/%s_strain_frequencies.pkl" % (config.analysis_directory,host,species))

    cluster_As = snp_alignment[:,:,0].T
    cluster_Ds = (snp_alignment.sum(axis=2)).T

    ## sort out indices with insufficient mean depth 
    good_inds = cluster_Ds.mean(axis=0) > 1
    good_inds = np.argwhere(good_inds == True).flatten()    

    snp_map = pd.read_pickle("%sclusters/%s/%s_snp_map.pkl" % (config.analysis_directory,host,species))
    snp_freqs = get_clusters_snv_trajectories(snp_map)
    
    strain_df = pd.DataFrame(columns=outp.keys())
    
    for K in outp.keys():
        strain_freq_est = outp[K]["centroid"][0]
        strain_df[K] = strain_freq_est
    
    strain_df.index = snp_samples
    
    dates = pd.read_pickle("%sPoyet_host_samples_timepoints.pkl" % config.metadata_directory)
    
    dates = dates[host]
    strain_df["Collection_Date"] = pd.Series(dates)
    strain_df = strain_df.sort_values("Collection_Date")
    
    strain_df["Order"] = range(strain_df.shape[0])
    sample_order = np.array(list(strain_df.loc[snp_samples]["Order"]))
    idx = np.empty_like(sample_order)
    idx[sample_order] = np.arange(len(sample_order))    
    strain_freqs = pd.DataFrame(index=strain_df.index,columns=outp.keys())
    strain_freqs[list(outp.keys())] = strain_df[outp.keys()]    
  #  strain_freqs = pd.DataFrame(index=strain_df.index,columns=outp.keys())

  #  strain_freqs[list(outp.keys())] = strain_df[outp.keys()]
    for key in snp_freqs.keys():
        snp_freqs[key] = snp_freqs[key].T[idx,:]
    
    return snp_freqs, strain_freqs,strain_df

def plot_freqs_snvs(snp_freqs,strain_freqs,strain_df,host,species,tpmin=0,tpmax=None,alpha=.1):
    
    import random
    import figure_utils
    
    colors = {}
    
    if host == "alien2" or host == "alien3" and strain_df.iloc[0]["Collection_Date"] == 0:
        tpmin = 1
    
    for key in list(strain_freqs.columns):
        r = random.random()
        b = random.random()
        g = random.random()
        colors[key] = (r, g, b)
    
    fig,ax = plt.subplots(1,figsize=(16,8))
    species_name = figure_utils.get_pretty_species_name(species)
    ax.set_title(species_name,size=25)
    
    for key in snp_freqs.keys():
        ax.plot(list(strain_df["Collection_Date"])[tpmin:tpmax],snp_freqs[key][:,:1000][tpmin:tpmax],color=colors[key],linewidth=.5,alpha=alpha);
        
        ax.plot(list(strain_df["Collection_Date"])[tpmin:tpmax],strain_freqs[key][tpmin:tpmax],color=colors[key],linewidth=2,linestyle="--");
        
        ax.scatter(list(strain_df["Collection_Date"])[tpmin:tpmax],strain_freqs[key][tpmin:tpmax],c=colors[key],edgecolors="k",s=100);
        
    if host == "alien":
        ax.axvspan(376, 380, alpha=0.5, color='red')
        ax.axvspan(630, 637, alpha=0.5, color='brown')
    
    fig.savefig("%sstrain_figures/%s/%s_strain_plot.png" % (config.analysis_directory,host,species))
    
def get_strain_total_freqs(species,host):
    
    anal_dir = config.analysis_directory
    
    oligo_species = config.oligo_species
    
    
    if species in oligo_species[host]:
        
        output_directory = "/u/scratch/r/rwolff/strainfinder_input/Poyet/%s" % host
        filename_prefix = "%s/%s" % (output_directory, species)
        snp_locations = pickle.load(open(filename_prefix+".strainfinder.locations.p",'rb'))
        snp_alignment = pd.read_pickle(filename_prefix+".strainfinder.p")
        snp_samples = pickle.load(open(filename_prefix+".strainfinder.samples.p",'rb'))
        snp_samples = [elem.decode("utf-8") for elem in snp_samples]

        dates = pd.read_pickle("metadata/Poyet_collection_dates.pkl")
        dates = pd.DataFrame(dates)
#dates = dates.loc[snp_samples]

        dates["Collection_Date"] = pd.to_datetime(dates.Collection_Date)
        outp=pd.read_pickle(f"{anal_dir}clusters/{host}/{species}_strain_frequencies.pkl")

        cluster_As = snp_alignment[:,:,0].T
        cluster_Ds = (snp_alignment.sum(axis=2)).T

        good_inds = cluster_Ds.mean(axis=0) > 1
        good_inds = np.argwhere(good_inds == True).flatten()

        snp_map = pd.read_pickle(f"~/diversity_ecology/analysis/clusters/{host}/{species}_snp_map.pkl")
        freqs = get_clusters_snv_trajectories(snp_map)
        
       # if freqs == {}:
       #     freqs[1] = 
        
        strain_df = pd.DataFrame(columns=outp.keys())
        
        for K in outp.keys():
            strain_freq_est = outp[K]["centroid"][0]
            strain_df[K] = strain_freq_est

        strain_df.index = snp_samples
        
        sra = pd.read_csv(f"{config.metadata_directory}/Poyet_SRA_report.txt",sep="\t",index_col=3)
        samps = sra["read_count"].loc[strain_df.index]
        
        ## filter out samples which have fewer than 1m total reads from further analysis
        samps = samps[samps > 1e6]
        strain_df = strain_df.loc[samps.index]

        strain_df["Collection_Date"] = dates["Collection_Date"]
        strain_df["Collection_Date"] = dates["Collection_Date"]
        strain_df = strain_df.sort_values("Collection_Date")
        strain_df["Date_Diffs"] = strain_df["Collection_Date"].diff().dt.days
        strain_df["Date_Diffs"] = strain_df["Date_Diffs"].replace(0.0,1)
        strain_df["Date_Diffs"][0] = 0.0
        strain_df["Order"] = range(strain_df.shape[0])
        sample_order = np.array(list(strain_df.loc[snp_samples]["Order"]))
        idx = np.empty_like(sample_order)
        idx[sample_order] = np.arange(len(sample_order))

        strain_freqs = pd.DataFrame(index=strain_df.index,columns=outp.keys())

        strain_freqs[list(outp.keys())] = strain_df[outp.keys()]
        out_strain_num = max(list(outp.keys())) + 1
        strain_freqs[out_strain_num] = 1 - strain_freqs.T.sum()
        times = np.array(list(strain_df["Date_Diffs"].cumsum()))

        rel_ab = pd.read_csv(f"{config.data_directory}species/relative_abundance.txt.bz2",sep="\t",index_col=0)
        spec_rel_ab = rel_ab.loc[species]
        spec_rel_ab = spec_rel_ab.loc[strain_df.index]
        strain_total_freqs = (strain_freqs.T*spec_rel_ab).T
        strain_total_freqs[["Date_Diffs","Order"]] = strain_df[["Date_Diffs","Order"]]
        strain_total_freqs["Date_Diffs"] = strain_total_freqs["Date_Diffs"].astype(int)
        
        ## again, safety check for this mislabeled sample
        if host == "ae":
            if "SRR9224093" in strain_total_freqs.index:
                strain_total_freqs = strain_total_freqs.drop("SRR9224093")
    
        return(strain_total_freqs)

    else:
        
        rel_ab = pd.read_csv(f"{config.data_directory}species/relative_abundance.txt.bz2",sep="\t",index_col=0)
        spec_rel_ab = rel_ab.loc[species]
        samples = config.Poyet_samples
        strain_total_freqs = pd.DataFrame(spec_rel_ab.loc[samples[host]])
        
        dates = pd.read_pickle("metadata/Poyet_collection_dates.pkl")
        dates = pd.DataFrame(dates)    

        dates["Collection_Date"] = pd.to_datetime(dates.Collection_Date)   
        
        strain_total_freqs["Collection_Date"] = dates["Collection_Date"]
        
        strain_total_freqs = strain_total_freqs.sort_values("Collection_Date")
        
        sra = pd.read_csv(f"{config.metadata_directory}/Poyet_SRA_report.txt",sep="\t",index_col=3)
        samps = sra["read_count"].loc[strain_total_freqs.index]
        ## filter by sample read depth
        samps = samps[samps > 1e6]
        strain_total_freqs = strain_total_freqs.loc[samps.index]
        
        strain_total_freqs["Date_Diffs"] = strain_total_freqs["Collection_Date"].diff().dt.days
        strain_total_freqs["Date_Diffs"] = strain_total_freqs["Date_Diffs"].replace(0.0,1)
        strain_total_freqs["Date_Diffs"][0] = 0.0
        strain_total_freqs["Order"] = range(strain_total_freqs.shape[0])
        
        strain_total_freqs = strain_total_freqs.rename(columns={species:1})
        strain_total_freqs["Date_Diffs"] = strain_total_freqs["Date_Diffs"].astype(int)
        
        del strain_total_freqs["Collection_Date"]
        
        return(strain_total_freqs)
         
def get_strain_total_freqs_Korpela(species,host):
    sra = pd.read_csv("~/diversity_ecology/scripts/metadata/Korpela_SRA_report.txt",sep="\t",
                index_col=5)
    output_directory = "/u/scratch/r/rwolff/strainfinder_input/Korpela/%s" % host
    filename_prefix = "%s/%s" % (output_directory, species)
    #snp_locations = pickle.load(open(filename_prefix+".strainfinder.locations.p",'rb'))
    snp_alignment = pd.read_pickle(filename_prefix+".strainfinder.p")
    snp_samples = pickle.load(open(filename_prefix+".strainfinder.samples.p",'rb'))
    snp_samples = [elem.decode("utf-8") for elem in snp_samples]
    sample_sra = sra.loc[[elem[1:] for elem in snp_samples]]
    dates = [d[2] for d in sample_sra["sample_alias"].str.split("-")]
    dates = {("c"+ind):int(sample_sra.loc[ind]["sample_alias"].split("-")[2]) for ind in sample_sra.index}
    outp=pd.read_pickle(f"~/diversity_ecology/analysis/clusters/{host}/{species}_strain_frequencies.pkl")

    cluster_As = snp_alignment[:,:,0].T
    cluster_Ds = (snp_alignment.sum(axis=2)).T
    
    good_inds = cluster_Ds.mean(axis=0) > 2
    good_inds = np.argwhere(good_inds == True).flatten()

    snp_map = pd.read_pickle(f"~/diversity_ecology/analysis/clusters/{host}/{species}_snp_map.pkl")
    freqs = get_clusters_snv_trajectories(snp_map)
    strain_df = pd.DataFrame(columns=outp.keys())
    for K in outp.keys():
        strain_freq_est = outp[K]["centroid"][0]
        strain_df[K] = strain_freq_est
    strain_df.index = snp_samples
    strain_df["Collection_Date"] = pd.Series(dates)
    strain_df = strain_df.sort_values("Collection_Date")
    strain_df["Order"] = range(strain_df.shape[0])
    sample_order = np.array(list(strain_df.loc[snp_samples]["Order"]))
    idx = np.empty_like(sample_order)
    idx[sample_order] = np.arange(len(sample_order))
    strain_freqs = pd.DataFrame(index=strain_df.index,columns=outp.keys())

    strain_freqs[list(outp.keys())] = strain_df[outp.keys()]
    out_strain_num = max(list(outp.keys())) + 1

    strain_freqs[out_strain_num] = 1 - strain_freqs.T.sum()
    rel_ab = pd.read_csv("/u/scratch/r/rwolff/merged_MIDAS_output/Korpela/species/relative_abundance.txt.bz2",sep="\t",index_col=0)
    spec_rel_ab = rel_ab.loc[species]
    spec_rel_ab = spec_rel_ab.loc[strain_df.index]
    strain_total_freqs = (strain_freqs.T*spec_rel_ab).T
    
    date_diffs = np.insert(np.diff(sorted(np.array(list(dates.values())))),0,0)
    
    strain_total_freqs["Date_Diffs"] = date_diffs
    strain_total_freqs["Order"] = range(0,strain_total_freqs.shape[0])
    
    cols = list(strain_total_freqs.columns)
    cols[:-2] = range(1,len(cols[:-2])+1)
    strain_total_freqs.columns = cols
   
                      
    return(strain_total_freqs)