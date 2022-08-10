def return_subject_sample_time_map():
    
    from collections import defaultdict

    HMP_directory = '/u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2'
    metadata_directory = '/u/home/r/rwolff/strain_stability_revisions/strainstability/metadata'

    sample_ids_file = open('%s/HMP1-2_samples.txt' % metadata_directory, 'r')
    metadata_file = open('%s/HMP1-2_metadata.txt' % metadata_directory, 'r')

    # Get all relevant sample IDs for this campaign
    all_run_accs = []
    for line in sample_ids_file:
            sample = line.strip()
            all_run_accs.append(sample)

    # Maps subject ID to dictionary of timepoint-run accession
    # (only for desired samples)
    subject_sample_time_map = defaultdict(dict)

    metadata_file.readline() # Ignore first line

    for line in metadata_file:
            subject_id, sample_id, run_acc, country, continent, visno = line.strip().split('\t')
            timept = visno
            if sample_id in all_run_accs:
                    subject_sample_time_map[subject_id][timept] = sample_id
            elif (sample_id + 'c') in all_run_accs:
                    subject_sample_time_map[subject_id][timept] = (sample_id + 'c')
    
    return subject_sample_time_map