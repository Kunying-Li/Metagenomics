## Installation

```
wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.2-py38-osx-conda.yml 
conda env create -n qiime2-amplicon-2024.2 --file qiime2-amplicon-2024.2-py38-osx-conda.yml  
conda activate qiime2-amplicon-2024.2
```

### prepare a data file for QIMME2

save it in sample_path.csv

```
sample-id,absolute-filepath,direction
sample-1,/Users/kunying/Downloads/Sam_data/1a_S1_L001_R1_001.fastq.gz,forward
sample-1,/Users/kunying/Downloads/Sam_data/1a_S1_L001_R2_001.fastq.gz,reverse
```

### import data to qza form

```
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path sample_path.csv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33
```

### visiualize sequence quality
```
qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization demux.qzv
```
upload the demux.qzv file to view.qiime2.org and click interactive quality plot to see the results

### run DATA2 to denoise data
```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 200 \
  --p-trunc-len-r 200 \
  --o-table table-dada2.qza \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-denoising-stats denoising-stats.qza

```
Here, you need to choose a value for the -p-trunc. Ask me or the AI about the details.

### taxonomy annotation
```
wget -O "silva-138-99-515-806-nb-classifier.qza" "https://data.qiime2.org/2024.2/common/silva-138-99-515-806-nb-classifier.qza"
```

```
qiime feature-classifier classify-sklearn \
  --i-classifier /Users/kunying/Downloads/Sam_data/silva-138-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs-dada2.qza \
  --o-classification taxonomy.qza
```
### taxonomy visualization
```
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```
### export data
```
qiime tools export \
  --input-path taxonomy.qza \
  --output-path result.csv

qiime tools export

```
### get the relative abundance
```
# Collapse to genus level (level 6 usually corresponds to genus level)
qiime taxa collapse \
  --i-table table-dada2.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table table-genus.qza

  
qiime feature-table relative-frequency \
  --i-table table-genus.qza \
  --o-relative-frequency-table table-genus-relative.qza

qiime tools export \
  --input-path table-genus-relative.qza \
  --output-path exported-table-genus-relative

biom convert -i exported-table-genus-relative/feature-table.biom -o table-genus-relative.tsv --to-tsv

```