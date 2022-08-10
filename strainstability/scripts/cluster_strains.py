import cluster_utils_exp as cu
import matplotlib.pyplot as plt
import cPickle
import pandas as pd
import state_utils
import config
import figure_utils as fu
import numpy
from itertools import cycle
from sklearn import metrics 
from math import log
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram,inconsistent,maxRstat
from scipy.spatial.distance import squareform
from parse_midas_data import parse_good_species_list
import sys

good_species = parse_good_species_list()

dates = config.dates
dates_diff = dates.diff()
spec_df = pd.read_csv("%sspecies/relative_abundance.txt.bz2" % config.data_directory,index_col=0,sep="\t")

strain_tot_df = pd.DataFrame(columns=spec_df.columns)

cluster_min_SNV_size = config.cluster_min_SNV_size
strainfinder_dir = config.strainfinder_directory 
dates = config.dates

min_coverage = 5

for species in good_species:
    base_dist_path = "/u/scratch/r/rwolff/distance_matrices/%s/%s/%s" % (config.cohort,config.host,species)

    sys.stderr.write("Processing %s \n\n" % species)

    snp_alignment = pd.read_pickle("%s/%s.strainfinder.p" %  (strainfinder_dir ,species))
    samples = pd.read_pickle("%s/%s.strainfinder.samples.p" % (strainfinder_dir ,species))
    samples = [s.decode("utf-8") for s in samples]

    snp_locations = pd.read_pickle("%s/%s.strainfinder.locations.p" % (strainfinder_dir,species))

    cluster_As = []
    cluster_Ds = []
    for snp_idx in range(0,snp_alignment.shape[1]):
        Ds = snp_alignment[:,snp_idx,:].sum(axis=1)
        As = snp_alignment[:,snp_idx,0]
        As = numpy.reshape(As, (1,len(As)))
        Ds = numpy.reshape(Ds, (1,len(Ds)))

        cluster_As.append(As[0])
        cluster_Ds.append(Ds[0])

    cluster_As = numpy.array(cluster_As)
    cluster_Ds = numpy.array(cluster_Ds)

    As = pd.DataFrame(cluster_As,columns=samples,index=snp_locations)
    Ds = pd.DataFrame(cluster_Ds,columns=samples,index=snp_locations)

    samples_sorted=list(dates.loc[samples].sort_values().index)

    As = As.loc[:,samples_sorted]
    Ds = Ds.loc[:,samples_sorted]

    F = As/Ds

    good_locs = pd.read_pickle("%s_locations.pkl" % base_dist_path)
    As = As.loc[good_locs.values]
    Ds = Ds.loc[good_locs.values]
    F = F.loc[good_locs.values]


    cluster_As = As.values
    cluster_Ds = Ds.values
    
    reader = pd.read_csv("%s.csv" % base_dist_path,chunksize=1000,index_col=0, iterator=True,low_memory=False)  
    df = pd.DataFrame([])
    for chunk in reader:
        df = df.append(chunk)
        
    reader1 = pd.read_csv("%s_1.csv" % base_dist_path,chunksize=1000,index_col=0, iterator=True,low_memory=False)  
    reader2 = pd.read_csv("%s_2.csv" % base_dist_path,chunksize=1000,index_col=0, iterator=True,low_memory=False)  

    df_1 = pd.DataFrame([])
    for chunk in reader1:
        df_1 = df_1.append(chunk)

    df_2 = pd.DataFrame([])
    for chunk in reader2:
        df_2 = df_2.append(chunk)
    
    distance_matrix = df.values
    distance_matrix_1 = df_1.values
    distance_matrix_2 = df_2.values
    
    Y = squareform(df.values)
    sys.stderr.write("SciPy hierarchical clustering...\n")
    Z =  linkage(Y, method='ward')
    sys.stderr.write("Done!\n")
    max_num_clusters = 4
    num_clusterss = numpy.arange(2,max_num_clusters+1)
    silhouette_scores = []
    for num_clusters in num_clusterss:

        nodes = fcluster(Z, int(num_clusters), criterion="maxclust")

        num_realized_clusters = len(set(nodes))

        if num_realized_clusters==1:
            S = 0
        else:
            S = metrics.silhouette_score(distance_matrix, nodes, metric = 'precomputed')

        silhouette_scores.append(S)

        print(num_clusters, num_realized_clusters, S)

    silhouette_scores = numpy.array(silhouette_scores)
    num_clusters = num_clusterss[silhouette_scores.argmax()]
    Smax = silhouette_scores.max()
    print(num_clusters, Smax)
    if Smax < 0:
        nodes = numpy.ones(distance_matrix.shape[0])
    else:
        nodes = fcluster(Z, num_clusters, criterion="maxclust")    
    
    cluster_snp_map = {}
    for snp_idx in range(0,len(nodes)):

        cluster_label = nodes[snp_idx]
        if cluster_label not in cluster_snp_map:
            cluster_snp_map[cluster_label] = []

        cluster_snp_map[cluster_label].append(snp_idx)

    snp_flip_map = {snp_idx: False for snp_idx in range(0,len(nodes))}

    cluster_fs_map = {}

    cluster_As_map = {}
    cluster_Ds_map = {}
    cluster_avg_fs_map = {}
    cluster_total_Ds_map = {}
    sss= []
    for cluster_label in cluster_snp_map.keys():

        anchor_idx = cluster_snp_map[cluster_label][0]

        cluster_As_map[cluster_label] = [cluster_As[anchor_idx,:]]
        cluster_Ds_map[cluster_label] = [cluster_Ds[anchor_idx,:]]

        if len(cluster_snp_map[cluster_label]) > 1:

            for snp_idx in cluster_snp_map[cluster_label][1:]:

                ## only keep SNVs where polarization is very clear
                if distance_matrix_1[anchor_idx,snp_idx]/distance_matrix[anchor_idx,snp_idx] < 3 and distance_matrix_2[anchor_idx,snp_idx]/distance_matrix[anchor_idx,snp_idx] < 3:
                    cluster_snp_map[cluster_label].remove(snp_idx)
                    continue

                target_As = cluster_As[snp_idx,:]
                target_Ds = cluster_Ds[snp_idx,:]

                if distance_matrix_2[anchor_idx,snp_idx] < distance_matrix_1[anchor_idx,snp_idx]:
                    # re-polarize
                    target_As = target_Ds-target_As
                    snp_flip_map[snp_idx] = not snp_flip_map[snp_idx]

                cluster_As_map[cluster_label].append(target_As)
                cluster_Ds_map[cluster_label].append(target_Ds)


        cluster_As_map[cluster_label] = numpy.array(cluster_As_map[cluster_label])
        cluster_Ds_map[cluster_label] = numpy.array(cluster_Ds_map[cluster_label])
        cluster_total_Ds_map[cluster_label] = cluster_Ds_map[cluster_label].sum(axis=0)
        cluster_avg_fs_map[cluster_label] = cluster_As_map[cluster_label].sum(axis=0)*1.0/(cluster_total_Ds_map[cluster_label]+(cluster_total_Ds_map[cluster_label]==0))

        # now polarize whole cluster if necessary
        if (cluster_avg_fs_map[cluster_label][0]+cluster_avg_fs_map[cluster_label][1])/2.0 > 0.5:
            cluster_avg_fs_map[cluster_label] = 1-cluster_avg_fs_map[cluster_label]
            cluster_As_map[cluster_label] = cluster_Ds_map[cluster_label] - cluster_As_map[cluster_label]
            for snp_idx in cluster_snp_map[cluster_label]:
                snp_flip_map[snp_idx] = not snp_flip_map[snp_idx]    

    # now write output
    cluster_map = {}
    for cluster_label in cluster_snp_map:

        cluster_map[cluster_label] = {}  
        cluster_map[cluster_label]['centroid'] = (cluster_avg_fs_map[cluster_label], cluster_total_Ds_map[cluster_label])
        cluster_map[cluster_label]['snps'] = {}
        for snp_idx in cluster_snp_map[cluster_label]:
            cluster_map[cluster_label]['snps'][snp_idx] = snp_flip_map[snp_idx]
        cluster_map[cluster_label]['snps'] = pd.Series(cluster_map[cluster_label]['snps'])

    for cluster_label in cluster_snp_map:
        cluster_map[cluster_label]["snp_freqs"]  = F.mask(Ds < min_coverage).iloc[cluster_map[cluster_label]['snps'].index]
        cluster_map[cluster_label]["snp_freqs"].index = cluster_map[cluster_label]['snps'].index

        to_flip=numpy.where(cluster_map[cluster_label]['snps'].values == False)
        not_flip=numpy.where(cluster_map[cluster_label]['snps'].values == True)
        to_flip_f = cluster_map[cluster_label]["snp_freqs"].iloc[to_flip]
        not_flip_f = cluster_map[cluster_label]["snp_freqs"].iloc[not_flip]
        df_flipped = pd.DataFrame(index=cluster_map[cluster_label]["snp_freqs"].index,columns=cluster_map[1]["snp_freqs"].columns)

        df_flipped.loc[to_flip_f.index] = to_flip_f
        df_flipped.loc[not_flip_f.index] = not_flip_f

        cluster_map[cluster_label]["snp_freqs"] = df_flipped

    ## remove clusters of less than 1000 SNVs
    for cluster_label in cluster_map:
        if len(cluster_map[cluster_label]["snps"]) < 1000:
            del cluster_snp_map[cluster_label]    
    
#     fig,ax = plt.subplots(figsize=(12,8))

#     for cluster_label in cluster_snp_map:

#         ax.plot(cluster_map[cluster_label]["snp_freqs"].sample(1000).T.values,color=cycol.next(),alpha=.003);
#         ax.plot(cluster_map[cluster_label]["centroid"][0],ls="--",lw=2,color=cycol.next(),label=cluster_label);
    
#     ax.plot((1-pd.DataFrame([cluster_map[c]["centroid"][0] for c in cluster_snp_map]).sum()).values,ls="--",lw=2,color=cycol.next())
#     ax.legend()
#     ax.grid(True)
    
    if len(cluster_snp_map.keys()) > 0:
        df_strain = pd.DataFrame(index=["%s_%s" % (species,c) for c in cluster_snp_map],columns=F.columns)

        for cluster_label in cluster_snp_map:
            df_strain.loc["%s_%s" % (species,cluster_label)] = cluster_map[cluster_label]["centroid"][0]
        
        max_ind = max([cluster_label for cluster_label in cluster_snp_map]) + 1
        df_strain.loc["%s_%s" % (species,max_ind)] = 1-df_strain.sum()  
    
    else:
        df_strain = pd.DataFrame(1,index=["%s_%s" % (species,1)],columns=F.columns)
        
    spec_strain_df = (df_strain*spec_df.loc[species])
    
    strain_tot_df = strain_tot_df.append(spec_strain_df)
    print(strain_tot_df)
      
strain_tot_df.to_csv("strains_df_%s.csv" % config.host)