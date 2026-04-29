# one-time: ensure bioconda is configured (skip if already done)
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels r
conda config --set channel_priority strict

# create env
mamba create -n dexseq_grch38 -y \
  python=3.10 \
  samtools parallel \
  htseq \
  r-base=4.3 r-essentials \
  bioconductor-dexseq bioconductor-deseq2 bioconductor-biocparallel \
  r-optparse r-data.table r-ggplot2 r-pheatmap

conda activate dexseq_grch38
which samtools
which dexseq_prepare_annotation.py
which dexseq_count.py
R --version