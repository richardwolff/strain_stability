import pandas as pd
import os

strain_stability_dir = "/u/home/r/rwolff/strain_stability_revisions/strainstability"
metadata_directory = os.path.expanduser("~/strain_stability_revisions/strainstability/metadata")

hosts = ["am","an","ae","ao"]

Poyet_fastq_dir = "/u/scratch/r/rwolff/Poyet_fastq_files/"
Poyet_data_directory = "/u/scratch/r/rwolff/Poyet_midas_output/merged_midas_output"
Poyet_dates_loc = "%s/metagenomic_scripts/metagenomic_data_files/Poyet_collection_dates.pkl" % strain_stability_dir
Poyet_dates = pd.read_pickle(Poyet_dates_loc)

suez_dir = "/u/project/ngarud/Garud_lab/suez_montassier/"
Suez_data_directory = "/u/project/ngarud/Garud_lab/suez_montassier/merged_midas_output/everything"
suez_snps_dir = "%smerged_midas_output/everything/snps" % suez_dir
Suez_metadata = pd.read_csv("%s/Suez_metadata.csv" % metadata_directory,index_col=0,sep='\t')
df_ind = Suez_metadata.groupby("Individual")
suez_host_samples = {i: df_ind.get_group(i).index for i in df_ind.groups.keys()}

cluster_distance_threshold_reads = 25
cluster_min_SNV_size = 500
