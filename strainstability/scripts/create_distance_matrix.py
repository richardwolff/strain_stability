import cluster_utils_exp as cu
import cPickle
import pandas as pd
import state_utils
import config
import sys
import numpy

species = sys.argv[1]
strainfinder_dir = config.strainfinder_directory 
dates = config.dates

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

sys.stderr.write("Data loaded \n")

if F.shape[0] > 10000:
    F = F.sample(10000)

sys.stderr.write("Data downsampled \n")

#F = F.loc[(F.mask(Ds < 5).T > 0).mean() > .25]

As = As.loc[F.index]
Ds = Ds.loc[F.index]

print(As.shape)

sys.stderr.write("Calculating distance matrix \n")

distance_matrix, distance_matrix_1, distance_matrix_2 = cu.calculate_distance_matrix(As.values,Ds.values)

df = pd.DataFrame(distance_matrix)
df_1 = pd.DataFrame(distance_matrix_1)
df_2 = pd.DataFrame(distance_matrix_2)

df.to_csv("/u/scratch/r/rwolff/distance_matrices/%s/%s/%s.csv" % (config.cohort,config.host,species))
df_1.to_csv("/u/scratch/r/rwolff/distance_matrices/%s/%s/%s_1.csv" % (config.cohort,config.host,species))
df_2.to_csv("/u/scratch/r/rwolff/distance_matrices/%s/%s/%s_2.csv" % (config.cohort,config.host,species))
pd.Series(F.index).to_pickle("/u/scratch/r/rwolff/distance_matrices/%s/%s/%s_locations.pkl" % (config.cohort,config.host,species))
