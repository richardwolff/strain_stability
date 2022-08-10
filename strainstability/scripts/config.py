###############################################################################
#
# Set up default source and output directories
#
###############################################################################
import os.path 
from math import log10
import pandas as pd
import state_utils

state_vars = state_utils.get_state_vars()
host = state_vars["DESIRED_HOST"]
db_type = state_vars["DB_TYPE"]
cohort = state_vars["DATA_COHORT"]

metadata_directory = os.path.expanduser("~/strain_stability_revisions/strainstability/metadata")

if cohort == "Poyet" and db_type == "isolate":
    data_directory = "/u/scratch/r/rwolff/Poyet_midas_output/merged_midas_output/%s/" % host
    midas_files_directory = "/u/scratch/r/rwolff/Poyet_midas_output/midas_files/%s/" % host
    midas_directory = "/u/project/ngarud/rwolff/midas_isolate_db_built/%s/" % host
    host_samples_file = "%s/host_samples/poyet_%s.txt" % (metadata_directory,host)
    f = open(host_samples_file, "r")
    host_samples = f.readlines()
    host_samples = [fl.strip() for fl in host_samples] 
    ## In order to align reads to isolates MIDAS-1.3.2 was used, which produces output files with slightly different
    ## configurationfrom MIDAS-mod. The following dictionary is used in post-processing steps to account for the differences in
    ## the file write convention 
    
    postproc_idx_dic = {"variant_type":-2, "gene_name":1,"ref_freq_file":"%ssnps/%s/snps_freq.txt.bz2"}
    MIDAS_ver = "1.3.2"
    
    strainfinder_directory = "/u/scratch/r/rwolff/strainfinder_input/Poyet/%s" % host
    
    all_dates = pd.read_csv("%s/host_dates/Poyet_dates.txt" % metadata_directory,index_col=0)
    dates = all_dates[all_dates["subject"] == host].sort_values("timepoint")["timepoint"]   
    
elif cohort == "Poyet" and db_type == "standard":
    data_directory = "/u/scratch/r/rwolff/Poyet_midas_output_spu/merged_midas_output/%s/" % host
    midas_files_directory = "/u/scratch/r/rwolff/Poyet_midas_output_spu/midas_files/%s/" % host
    midas_directory = "/u/project/ngarud/Garud_lab/midas_db_v1.2/"  
    host_samples_file = "%s/host_samples/poyet_%s.txt" % (metadata_directory,host)
    f = open(host_samples_file, "r")
    host_samples = f.readlines()
    host_samples = [fl.strip() for fl in host_samples] 
    
    postproc_idx_dic = {"variant_type":5, "gene_name":6,"ref_freq_file":"%ssnps/%s/snps_ref_freq.txt.bz2"}
    MIDAS_ver = "mod"
    
    strainfinder_directory = "/u/scratch/r/rwolff/strainfinder_input/Poyet/%s" % host
    
    all_dates = pd.read_csv("%s/host_dates/Poyet_dates.txt" % metadata_directory,index_col=0)
    dates = all_dates[all_dates["subject"] == host].sort_values("timepoint")["timepoint"]
    
elif cohort == "Suez":
    data_directory = "/u/project/ngarud/Garud_lab/suez_montassier/merged_midas_output/host_merged/%s/" % host
    midas_files_directory = "/u/project/ngarud/Garud_lab/suez_montassier/midas_files/everything/"
    midas_directory = "/u/project/ngarud/Garud_lab/midas_db_v1.2/"
    host_samples_file = "%s/host_samples/suez_%s.txt" % (metadata_directory,host)
    f = open(host_samples_file, "r")
    host_samples = f.readlines()   
    host_samples = [fl.strip() for fl in host_samples] 

    postproc_idx_dic = {"variant_type":5, "gene_name":6,"ref_freq_file":"%ssnps/%s/snps_ref_freq.txt.bz2"}
    MIDAS_ver = "mod"

    strainfinder_directory = "/u/scratch/r/rwolff/strainfinder_input/Suez/%s" % host
    
elif cohort == "Roodgar":
    data_directory = "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Roodgar_Good_data/highres_microbiome_timecourse_data_022119/"
    midas_files_directory = None
    midas_directory = None
    host_samples_file = "%s/host_samples/Roodgar.txt" % (metadata_directory)
    f = open(host_samples_file, "r")
    host_samples = f.readlines()   
    host_samples = [fl.strip() for fl in host_samples] 

    postproc_idx_dic = {"variant_type":5, "gene_name":6,"ref_freq_file":"%ssnps/%s/snps_ref_freq.txt.bz2"}
    
    strainfinder_directory = "/u/scratch/r/rwolff/strainfinder_input/Roodgar"
    host = None
    
    dates = pd.read_csv("%s/host_dates/Roodgar_dates.txt" % metadata_directory,index_col=1)["timepoint"]
    
elif cohort == "Korpela":    
    data_directory = "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Korpela/merged_midas_output/host_merged/%s/" % host
    midas_files_directory =  "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Korpela/midas_output/"
    midas_directory = "/u/project/ngarud/Garud_lab/midas_db_v1.2/"
    host_samples_file = "%s/host_samples/Korpela_%s.txt" % (metadata_directory,host)
    f = open(host_samples_file, "r")
    host_samples = f.readlines()   
    host_samples = [fl.strip() for fl in host_samples] 

    MIDAS_ver = "mod"
    postproc_idx_dic = {"variant_type":5, "gene_name":6,"ref_freq_file":"%ssnps/%s/snps_ref_freq.txt.bz2"}

    strainfinder_directory = "/u/scratch/r/rwolff/strainfinder_input/Korpela/%s" % host
    
    all_dates = pd.read_csv("%s/host_dates/Korpela_dates.txt" % metadata_directory,index_col=0,sep="\t")
    dates = all_dates[all_dates["subject"] == host].sort_values("timepoint")["timepoint"]

    
elif cohort == "Schirmer":    
    data_directory = "/u/scratch/r/rwolff/Schirmer_midas_output/merged_midas_output/%s/" % host
    midas_files_directory = "/u/scratch/r/rwolff/Schirmer_midas_output/midas_files/%s/" % host 
    midas_directory = "/u/project/ngarud/Garud_lab/midas_db_v1.2/"

    host_samples_file = "%s/host_samples/Schirmer_%s.txt" % (metadata_directory,host)
    f = open(host_samples_file, "r")
    host_samples = f.readlines()   
    host_samples = [fl.strip() for fl in host_samples] 

    MIDAS_ver = "mod"
    postproc_idx_dic = {"variant_type":5, "gene_name":6,"ref_freq_file":"%ssnps/%s/snps_ref_freq.txt.bz2"}

    strainfinder_directory = "/u/scratch/r/rwolff/strainfinder_input/Schirmer/%s" % host

## TKTK
    all_dates = pd.read_csv("%s/host_dates/Schirmer_dates.txt" % metadata_directory,index_col=0)
    dates = all_dates[all_dates["subject"] == host].sort_values("timepoint")["timepoint"]

elif cohort == "HMP":    
    data_directory = "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/data/"
    midas_files_directory = ""
    midas_directory = ""
    host_samples_file = ""
    host_samples = ""
    ## In order to align reads to isolates MIDAS-1.3.2 was used, which produces output files with slightly different
    ## configurationfrom MIDAS-mod. The following dictionary is used in post-processing steps to account for the differences in
    ## the file write convention 
    
    postproc_idx_dic = ""
    MIDAS_ver = ""
    
    strainfinder_directory = ""
    
    all_dates = ""
    dates = "" 
    
    import HMP_time_pair_species_list
    subject_sample_time_map = HMP_time_pair_species_list.return_subject_sample_time_map()
    
analysis_directory = os.path.expanduser("~/strain_stability_revisions/strainstability/analysis")
scripts_directory = os.path.expanduser("~/strain_stability_revisions/strainstability/scripts/")
patric_directory = os.path.expanduser("/u/project/ngarud/Garud_lab/software/PATRIC/")

# We use this one to debug because it was the first one we looked at
debug_species_name = 'Bacteroides_uniformis_57318'

## what is the minimum depth a site must have to be included in pi calculations? 
min_depth_pi_site = 5
min_prev_pi_site = .2

good_species_min_coverage = 2
good_species_min_prevalence = 2

min_median_coverage = 20

consensus_lower_threshold = 0.2
consensus_upper_threshold = 0.8
fixation_min_change = consensus_upper_threshold-consensus_lower_threshold
fixation_log10_depth_ratio_threshold = log10(3)

gainloss_max_absent_copynum = 0.05
gainloss_min_normal_copynum = 0.6
gainloss_max_normal_copynum = 1.2

core_genome_min_copynum = 0.3
core_genome_max_copynum = 3 #
core_genome_min_prevalence = 0.9
shared_genome_min_copynum = 3

# Default parameters for pipe snps
# (Initial filtering for snps, done during postprocessing)
pipe_snps_min_samples=4
pipe_snps_min_nonzero_median_coverage=5
pipe_snps_lower_depth_factor=0.3
pipe_snps_upper_depth_factor=3

parse_snps_min_freq = 0.05

max_d = 3.5
#cluster_min_coverage = 5
cluster_min_coverage = 10
cluster_min_SNV_size = 1000


