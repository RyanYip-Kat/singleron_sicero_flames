#/home/ye/Work/BioAligment/Nanopore/FLAMES/src/bin/match_cell_barcode /home/ye/Work/BioAligment/Nanopore/singleron/example_data/ example.log example_singleron.fastq.gz  /home/ye/Work/BioAligment/Nanopore/singleron/example_sicelore_barcode.tsv 3 45

genomeFa="/home/ye/Data/10X/VDJ/ref/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa"
outdir="/home/ye/Work/BioAligment/Nanopore/FLAMES/example_singleron"
gtf="/home/ye/Data/10X/VDJ/ref/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf"
#gtf="GRCh38.gff3"
python/sc_long_pipeline.py -a $gtf -i example_singleron.fastq.gz --minimap2_dir /home/ye/anaconda3/envs/FLAMES/bin/ --genomefa $genomeFa --outdir $outdir --downsample_ratio 0.25
