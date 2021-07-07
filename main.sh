data_dir="/home/ye/Work/BioAligment/Nanopore/singleron/scONT/data/"
genomeFa="/home/ye/Data/10X/VDJ/ref/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa"
gtf="/home/ye/Data/10X/VDJ/ref/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf"
outdir="/home/ye/Work/BioAligment/Nanopore/FLAMES/singleron_result"
whitelist_barcode="/home/ye/Work/BioAligment/Nanopore/FLAMES/barcodes.tsv"

/home/ye/Work/BioAligment/Nanopore/FLAMES/src/bin/match_cell_barcode $data_dir $outdir/run.log $outdir/match_singleron.fastq.gz $whitelist_barcode 3 45

python/sc_long_pipeline.py -a $gtf -i $outdir/match_singleron.fastq.gz  --minimap2_dir /home/ye/anaconda3/envs/FLAMES/bin/ --genomefa $genomeFa --outdir $outdir 
