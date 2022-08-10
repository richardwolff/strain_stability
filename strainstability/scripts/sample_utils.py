import numpy as np
import config, parse_midas_data
from collections import defaultdict

# ===========================================================================
# This module (parse_metadata) contains the following utilities:
# 
# 	get_sample_names, calculate_qp_samples
# 	calculate_ordered_same_subject_pairs, calculate_ordered_diff_subject_pairs
# 	calculate_ordered_same_sample_pairs, calculate_ordered_subject_pairs
# 	calculate_sample_idx_map
# 	apply_sample_index_map_to_indices, calculate_samples_in_different_subjects
# 	calculate_subject_pairs, calculate_sample_subject_matrix
# 
# 	parse_sample_metadata_map, extract_sample_metadata_map
# 	parse_sample_subject_map, parse_sample_country_map, parse_sample_continent_map
# 	parse_subject_sample_map, parse_sample_order_map
# 
# 	parse_merged_sample_names, calculate_unique_samples
# ===========================================================================

# Generic helper functions

flatten = lambda l: [item for sublist in l for item in sublist]

# ===========================================================================
# FUNCTIONS FOR DEALING WITH SAMPLE LISTS
# ===========================================================================

# ===========================================================================
# get_sample_names: returns list of samples in specified host
# 
# host: one of 'all', 'ae', 'ao', 'am', 'an'
# ===========================================================================

def get_sample_names(desired_host = 'all'):
	
	sample_dict = defaultdict(set) # Keys are host IDs (ae, ao, am, an)
	all_samples = []
	
	# Poyet
	samples_fpath = "%s/Poyet_run_accessions_actual_alias.txt" % (config.metadata_directory)
	
	with open(samples_fpath, 'r') as file:
		header = file.readline()
		for line in file:
			sample, alias = line.strip().split('\t')
			host, tp = alias.split('-')
			sample_dict[host].add(sample)
			all_samples.append(sample)
	
	return (all_samples if desired_host == 'all' else sample_dict[desired_host])

# ===========================================================================
# calculate_qp_samples: returns QP status dictionary given samples, species
# 'qp', 'non-qp' [high coverage], 'low-coverage'
# ===========================================================================

def calculate_qp_samples(all_samples, species_name):
	
	import diversity_utils
	
	# list of samples that meet coverage criteria for this species
	highcoverage_samples = set(diversity_utils.calculate_highcoverage_samples(species_name))
	
	# list of samples that meet QP criteria for this species
	haploid_samples = set(diversity_utils.calculate_haploid_samples(species_name))
	
	qp_statuses = ['qp','non-qp','low-coverage']
	qp_sample_sets = {status: set() for status in qp_statuses}
	
	for sample in all_samples:
		if sample in haploid_samples:	# QP (must be high coverage)
			qp_sample_sets['qp'].add(sample)
		elif sample in highcoverage_samples:	# Non-QP (high coverage)
			qp_sample_sets['non-qp'].add(sample)
		else:	# Low coverage
			qp_sample_sets['low-coverage'].add(sample)
	
	return qp_sample_sets

# ===========================================================================
# load_qp_samples: returns QP status dictionary given samples, species
# 'qp', 'non-qp' [high coverage], 'low-coverage' (uses pickle)
# ===========================================================================

def load_qp_samples(desired_samples, species_name, force_repickle = False):
	
	import pickle, os.path
	pickle_fn = "%s/pickles/qp_samples/%s_qp_sample_dict.pkl" % (config.data_directory, species_name)
	
	if force_repickle or not os.path.isfile(pickle_fn):
		all_samples = get_sample_names('all')
		qp_sample_sets = calculate_qp_samples(all_samples, species_name)
		pickle.dump(qp_sample_sets, open(pickle_fn, 'wb'))
		return qp_sample_sets
	else:
		qp_sample_sets = pickle.load(open(pickle_fn, 'rb'))
		for cat in qp_sample_sets: # qp, non-qp, low-coverage
			old_sample_list = list(qp_sample_sets[cat])
			for sample in old_sample_list:
				if sample not in desired_samples:
					qp_sample_sets[cat].remove(sample)
		return qp_sample_sets

# ===========================================================================
# FUNCTIONS FOR DEALING WITH SAMPLE PAIRS
# ===========================================================================

# ===========================================================================
# calculate_ordered_same_subject_pairs: computes same-subject sample pairs;
# specify if timepoints should be nonconsecutive, consecutive, or longest
#
# Returns same_subject_idxs, tuple with idx1 (lower order) and idx2 (higher)
# ===========================================================================

def calculate_ordered_same_subject_pairs(sample_order_map, sample_list=[], within_host_type='consecutive'):
	
	idx_lower, idx_upper = [], []
	
	# Reconstruct "timeseries" for each subject
	subject_order_idx_map = defaultdict(dict)
	for i in range(0,len(sample_list)):
		subject, order = sample_order_map[sample_list[i]]		 
		subject_order_idx_map[subject][order] = i
	
	# Create index pairs within subjects
	for subject in subject_order_idx_map:
		sorted_orders = list(sorted(subject_order_idx_map[subject].keys()))
		
		if len(sorted_orders) < 2:	# no pairs can be formed
			continue
		
		if within_host_type=='longest':
			idx_lower.append(subject_order_idx_map[subject][sorted_orders[0]])
			idx_upper.append(subject_order_idx_map[subject][sorted_orders[-1]])
		elif within_host_type=='consecutive':
			for order_idx in range(1,len(sorted_orders)):
				idx_lower.append(subject_order_idx_map[subject][sorted_orders[order_idx-1]])
				idx_upper.append(subject_order_idx_map[subject][sorted_orders[order_idx]])
		elif within_host_type=='nonconsecutive':
			for order_idx_i in range(0,len(sorted_orders)):
				for order_idx_j in range(order_idx_i+1,len(sorted_orders)):
					idx_lower.append(subject_order_idx_map[subject][sorted_orders[order_idx_i]])
					idx_upper.append(subject_order_idx_map[subject][sorted_orders[order_idx_j]])
	
	same_subject_idxs = (np.array(idx_lower), np.array(idx_upper))
	
	return same_subject_idxs

# ===========================================================================
# calculate_ordered_diff_subject_pairs: computes diff-subject sample pairs;
# only one sample considered for each subject; specify whether timepoint
# should be first or last # TODO?
#
# Returns diff_subject_idxs, tuple with idx1 and idx2
# ===========================================================================

def calculate_ordered_diff_subject_pairs(sample_order_map, sample_list=[], diff_host_type='first'):
	
	# Reconstruct "timeseries" for each subject
	subject_order_idx_map = defaultdict(dict)
	for i in xrange(0,len(sample_list)):
		subject, order = sample_order_map[sample_list[i]]		 
		subject_order_idx_map[subject][order] = i
	
	sorted_subjects = sorted(subject_order_idx_map.keys()) # all subjects
	op = max if diff_host_type == 'last' else min
	idx_lower, idx_upper = [], []
	
	for subject_i_idx in xrange(0,len(sorted_subjects)):
		subject_i = sorted_subjects[subject_i_idx]
		# Get earliest timepoint for this subject
		earliest_order_i = op(subject_order_idx_map[subject_i].keys())
		i = subject_order_idx_map[subject_i][earliest_order_i]
		
		# Pair up with earliest timepoint from all other subjects
		for subject_j_idx in xrange(subject_i_idx+1,len(sorted_subjects)):
			subject_j = sorted_subjects[subject_j_idx]
			earliest_order_j = op(subject_order_idx_map[subject_j].keys())
			j = subject_order_idx_map[subject_j][earliest_order_j]
			
			idx_lower.append(i)
			idx_upper.append(j)
	
	diff_subject_idxs = (np.array(idx_lower,dtype=np.int32), np.array(idx_upper,dtype=np.int32))
	return diff_subject_idxs

# ===========================================================================
# calculate_ordered_same_sample_pairs: computes same sample "pair" indices
# Assumes no duplicate samples in sample_list
# ===========================================================================

def calculate_ordered_same_sample_pairs(sample_order_map, sample_list=[]):
	idxs = np.arange(0,len(sample_list))
	return (np.array(idxs,dtype=np.int32), np.array(idxs,dtype=np.int32))

# ===========================================================================
# calculate_ordered_subject_pairs: wrapper function for combining 
# calculate_ordered_same_sample_pairs, calculate_ordered_same_subject_pairs
# and calculate_ordered_diff_subject_pairs
# ===========================================================================

def calculate_ordered_subject_pairs(sample_order_map, sample_list=[], within_host_type='consecutive', diff_host_type='first'):
	
	same_sample_idxs = calculate_ordered_same_sample_pairs(sample_order_map, sample_list)
	same_subject_idxs = calculate_ordered_same_subject_pairs(sample_order_map, sample_list, within_host_type)
	diff_subject_idxs = calculate_ordered_diff_subject_pairs(sample_order_map, sample_list, diff_host_type)
	return same_sample_idxs, same_subject_idxs, diff_subject_idxs

# ===========================================================================
# calculate_sample_idx_map: creates a map of indexes from one list of samples
# (sample_list_from) to another list of samples (sample_list_to).
# The from list must be a strict subset of the to list. 
# ===========================================================================

def calculate_sample_idx_map(sample_list_from, sample_list_to):
	
	sample_list_to = list(sample_list_to)
	sample_map = {}
	for i in xrange(0,len(sample_list_from)):
		sample_map[i] = sample_list_to.index(sample_list_from[i])
	
	return sample_map

def apply_sample_index_map_to_indices(sample_idx_map, idxs):
	new_idxs = (np.array([sample_idx_map[i] for i in idxs[0]]), np.array([sample_idx_map[i] for i in idxs[1]]))
	return new_idxs

# ===========================================================================
# calculate_samples_in_different_subjects: returns boolean array indicating
# whether each sample in sample_list has different subject as focal_sample
# ===========================================================================

def calculate_samples_in_different_subjects(sample_subject_map, sample_list, focal_sample):
	focal_subject = sample_subject_map[focal_sample]
	subjects = np.array([sample_subject_map[s] for s in sample_list])
	return (subjects != focal_subject)

# ===========================================================================
# calculate_subject_pairs: calculates which samples belong to different 
# subjects, which belong to different timepoints in same subject, and 
# which are the same timepoint.
#
# Returns same_sample_idxs, same_subject_idxs, diff_subject_idxs, 
# each of which is a tuple with lists idx1 and idx2.
# All pairs are included only once (order doesn't matter).
# ===========================================================================

def calculate_subject_pairs(sample_subject_map, sample_list = None):
	
	if sample_list is None:
		sample_list = sample_subject_map.keys()
	
	same_sample_idx_lower, same_sample_idx_upper = [], []
	same_subject_idx_lower, same_subject_idx_upper = [], []
	diff_subject_idx_lower, diff_subject_idx_upper = [], []
	
	for i in xrange(0,len(sample_list)):
		sample_i = sample_list[i]
		same_sample_idx_lower.append(i)
		same_sample_idx_upper.append(i)
		for j in xrange(0,i):
			sample_j = sample_list[j]
			if sample_subject_map[sample_i] == sample_subject_map[sample_j]:
				same_subject_idx_lower.append(i)
				same_subject_idx_upper.append(j)
			else: 
				diff_subject_idx_lower.append(i)
				diff_subject_idx_upper.append(j)
		
	same_sample_idxs = (np.array(same_sample_idx_lower,dtype=np.int32), np.array(same_sample_idx_upper,dtype=np.int32))
	same_subject_idxs = (np.array(same_subject_idx_lower,dtype=np.int32), np.array(same_subject_idx_upper,dtype=np.int32))	
	diff_subject_idxs = (np.array(diff_subject_idx_lower,dtype=np.int32), np.array(diff_subject_idx_upper,dtype=np.int32))
	
	return same_sample_idxs, same_subject_idxs, diff_subject_idxs

# ===========================================================================
# calculate_sample_subject_matrix: matrix, rows are subjects, columns are hosts 
# A_ih = 1 if sample i is in host h 
# ===========================================================================

def calculate_sample_subject_matrix(samples):
	
	sample_idx_map = {samples[i]:i for i in xrange(0,len(samples))}
	
	subject_sample_map = parse_subject_sample_map()
	subjects = subject_sample_map.keys()
	
	sample_subject_matrix = np.zeros((len(samples),len(subjects)),dtype=np.bool)
	
	for subject_idx in xrange(0,len(subjects)):
		for sample in subject_sample_map[subjects[subject_idx]]:
			if sample in sample_idx_map:
				sample_subject_matrix[sample_idx_map[sample], subject_idx] = True
	
	return sample_subject_matrix, subjects

# ===========================================================================
# FUNCTIONS FOR SAMPLE-METADATA MAPS
# ===========================================================================

# ====================================================================================
# parse_sample_read_count_map
# 
# Loads HMP read counts (# lines in fastq files / 4)
# Returns map: sample -> read_count
# ====================================================================================

def parse_sample_read_count_map(): 
	
	from config import metadata_directory
	sample_rc_map = {}
	
	# First load HMP metadata
	with open(metadata_directory + "Poyet_sample_read_counts.txt", 'r') as file:	
		for line in file:
			sample_id, read_count = line.strip().split('\t')
			sample_rc_map[sample_id] = int(read_count)
	
	return sample_rc_map

# ====================================================================================
# parse_sample_day_map
# 
# Loads collection date and uses that to compute number of days since first sample
# for each longitudinally sampled host
# ====================================================================================

def parse_sample_day_map():
	import csv
	from datetime import datetime
	
	days_between = lambda dt1, dt2: abs((dt1-dt2).days)
	
	reader = csv.reader(open('%s/Poyet_metadata_full.csv' % (config.metadata_directory), 'r'))
	rows = [row for row in reader]
	header = rows[0]

	sample_date_map = {}

	for row in rows[1:]:
		sample_acc = row[0]
		date = row[7]
		sample_date_map[sample_acc] = date
	
	subject_sample_map = parse_subject_sample_map()
	sample_day_map = {}
	
	for subject in subject_sample_map:		
		min_dt = min([datetime.strptime(sample_date_map[sample], "%d-%b-%Y") for sample in subject_sample_map[subject]])
		for sample in subject_sample_map[subject]:
			dt = datetime.strptime(sample_date_map[sample], "%d-%b-%Y")
			days = days_between(min_dt, dt)
			sample_day_map[sample] = days
	
	return sample_day_map

# ====================================================================================
# parse_sample_metadata_map
# 
# Loads metadata for Poyet samples
# Returns map:
# sample -> (subject_id, run_accession, temporal_order)
# ====================================================================================

def parse_sample_metadata_map(): 
	
	sample_metadata_map = {}
	
	# Load Poyet metadata
	fpath = '%s/Poyet_run_accessions_actual_alias.txt' % (config.metadata_directory)
	
	with open(fpath, 'r') as file:
		header = file.readline()
		for line in file:
			run_accession, alias = line.strip().split('\t')
			subject_id, suffix = alias.split('-')
			order = int(suffix.split('_')[0])
			sample_metadata_map[run_accession] = (subject_id, run_accession, order)
	
	return sample_metadata_map

# ====================================================================================
# extract_sample_metadata_map
# 
# Using passed in sample-metadata map, extracts information for only one column
# (options: subject_id, country, continent or order) and returns new map
# Loads the default full sample-metadata map if nothing passed in
# ====================================================================================

def extract_sample_metadata_map(field, sample_metadata_map = None):
	
	field_dict = {"subject_id": 0, "order": 2}
	field_idx = field_dict[field] if (field in field_dict) else 0 # Defaults to subject_id
	
	if sample_metadata_map is None:
		sample_metadata_map = parse_sample_metadata_map() # Load it
	
	extracted_sample_metadata_map = {}
	for sample in sample_metadata_map:
		extracted_sample_metadata_map[sample] = sample_metadata_map[sample][field_idx]
	
	return extracted_sample_metadata_map

# ====================================================================================
# parse_sample_subject_map, parse_sample_country_map, parse_sample_continent_map
#
# Convenience functions for extract_sample_metadata_map
# ====================================================================================

def parse_sample_subject_map(sample_metadata_map = None): 
	return extract_sample_metadata_map("subject_id", sample_metadata_map)

# ====================================================================================
# parse_subject_sample_map
# 
# Returns map: subject -> set of accession IDs
# ====================================================================================

def parse_subject_sample_map(sample_metadata_map = None): 
		
	if sample_metadata_map is None:
		sample_metadata_map = parse_sample_metadata_map() # Load it
	
	subject_sample_map = defaultdict(set)
	for sample in sample_metadata_map:
		subject_id, _, order = sample_metadata_map[sample]		
		subject_sample_map[subject_id].add(sample)
	
	return subject_sample_map

# ====================================================================================
# parse_sample_order_map
# 
# Returns map from sample -> (subject_id, temporal_order)
# ====================================================================================

def parse_sample_order_map(sample_metadata_map = None): 
	
	if sample_metadata_map is None:
		sample_metadata_map = parse_sample_metadata_map() # Load it
	
	sample_order_map = {}
	for sample in sample_metadata_map:
		subject_id, _, order = sample_metadata_map[sample]
		sample_order_map[sample] = (subject_id, order)
	
	return sample_order_map

# ===========================================================================
# Simply removes 'c' suffix for any merged samples
# (Unnecessary but leave here for compatibility)
# ===========================================================================

def parse_merged_sample_names(items):
		samples = []
		for item in items:
				sample = item.strip()
				if sample.endswith('c'):
						sample = sample[:-1]
				samples.append(sample)
		
		samples = np.array(samples)
		return samples

# ===========================================================================
# Prunes sample list to remove multiple timepoints from same subject
# Considers mothers and infants as different subjects
# Returns len(sample_list) boolean array with element=False if sample was pruned	
# ===========================================================================

def calculate_unique_samples(subject_sample_map, sample_list=[]):
	
	# If nothing passed in, default to all samples
	if len(sample_list)==0:
			sample_list = sorted(get_sample_names('all'))
	
	# Get sample subject map
	sample_subject_map = parse_sample_subject_map()
	
	subject_idx_map = {}
	
	for i in range(0,len(sample_list)):
		sample = sample_list[i]
		subject = sample_subject_map[sample]
		if subject not in subject_idx_map:
			subject_idx_map[subject] = i # Take the first sample of that subject
	
	unique_idxs = np.zeros(len(sample_list),dtype=np.bool_)
	for i in subject_idx_map.values():
			unique_idxs[i]=True
	
	return unique_idxs


