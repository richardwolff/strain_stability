export PATH=$PATH:/u/project/ngarud/rwolff/MIDAS/scripts
export PYTHONPATH=/u/project/ngarud/rwolff/MIDAS
export PATH=$PATH:/u/project/ngarud/rwolff/hmmer-3.3.2
export PATH=$PATH:/u/project/ngarud/rwolff/vsearch-2.21.1

host=$1
base_dir=/u/project/ngarud/rwolff/midas_isolate_db
indir=$base_dir/$host
mapfile=$base_dir/$host.mapfile
outdir=/u/project/ngarud/rwolff/midas_isolate_db_built/$host

build_midas_db.py $indir $mapfile $outdir

## note: uncommented line to build genome.features