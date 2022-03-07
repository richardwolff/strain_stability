. /u/local/Modules/default/init/modules.sh
module load ncbi

readarray accessions < /u/home/NIH1/Olm_NIH1_run_accessions.txt

i=$1
accession=${accessions[$i]}
echo $accession

fastq_dir=/u/scratch/r/rwolff/Poyet_fastq_files

SRA_fpath=/u/scratch/r/rwolff/Poyet_sra_files/${accession}

fastq-dump-orig.2.10.1 -O $fastq_dir $SRA_fpath --gzip --split-files
