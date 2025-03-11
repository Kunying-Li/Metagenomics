
# Metagenomics Data Processing Protocol

This document provides an example of how to process metagenomics raw data and build metagenomics assemblies. The protocol is divided into four main parts, each requiring the use of specific software tools.

## Software Installation Guidelines

For each part of the protocol, you will need to utilize multiple software tools. It is typical to download the source code from GitHub, compile it in your personal folder, and call them with the installation directory, especially when using the MIT Engaging cluster. Using pip or conda is not recommended on the MIT Engaging cluster because you would need to install them in your personal folder and manually add the environmental path.

## Protocol Steps

1. **Initial Filtering**: The first step involves filtering the raw data to remove low-quality reads. This step also includes calculating quality scores to assess the data.
   
2. **Sequence Alignment**: In this example, our samples consist of shrimp DNA and bacteria DNA. The goal is to remove all DNA sequences that match the reference DNA from *Litopenaeus vannamei*.
   
3. **Building Contigs**: The third step focuses on constructing contigs from single reads.
   
4. **Assemble genomes**: Fourth, we will build metagenome assemblies from the contigs using three different software, and do consolidation to get the most-likely metagenome data.
5. **Post-processing**: Finally, we need to do evaluation and other additional stuffs.

## About the MIT Engaging Cluster

Currently, we utilize the shared resources from the Chisholm lab. You will need to request permission to use their nodes. While it is possible to run processes on the login nodes, programs that are found to be computationally expensive will be terminated without notice.

This is the command to see how many nodes are used. 

```
sinfo -N -p sched_mit_chisholm -l
squeue -p sched_mit_chisholm -o "%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R"
```

This is the code to call just one node. You will find your server changed to Nodexxx after running it in the command line, you are allocated all the cpus and memory within one node for a day

```
salloc -p sched_mit_chisholm -n 20 --mem=250GB -t 1-00:00:00 
```

To submit a job, you need to create a .bash file with the following header. Write your own code after the header.

```
#!/bin/bash
#SBATCH -N 1 ## 
#SBATCH -n 20 ## cores, chisholm have 10 nodes each with 250 GB ram and 20 cpus/cores
#SBATCH -t 10-00:00:00 ## time limit, max in Chisholm is 90 days
#SBATCH -p sched_mit_chisholm ## tell it which scheduler
#SBATCH -J skewer_SAMPLE ## job name, unique for each one
#SBATCH --mem=250G ## if mapping, go a little bit above the database size
#SBATCH --job-name=example_job
#SBATCH --output=/path/to/logfile_%j.out
#SBATCH --error=/path/to/logfile_%j.err
```
Useful commands

```
sbatch test.bash ## Submission of the job
scancel #job_number ## Cancel the job
sacct ## check job status
```

Sometimes if you forgot to specify the name of the log file and miss up a bunch of log file and output file, it might be helpful to sort all the files by modification time to find the output file corresponding with its log file.

```
find . -type f -printf '%T+ %p\n' | sort -r
```

## Initial filtering

### Install trimmomatic
You will need to load jdk module before compiling it. Check the available module version if my code doesn't work. 

```
wget https://github.com/timflutre/trimmomatic/archive/refs/heads/master.zip
module add jdk/18.0.1.1
salloc -p sched_mit_chisholm -n 20 --mem=250GB -t 1-00:00:00 ##compiling this module is a little computational expensive, call a node to do that
make ##cd into the main folder
```
### Example code and parameter choosing
The following code is given the first two file as the input, and output trimmed sequence. The parameter is a set of common settings and you can check with GPT to see the details.

```
java -jar /home/kunyin71/software/trimmomatic/classes/trimmomatic.jar PE \
/nobackup1/kunyin71/Illumina_DNA_Reads/raw/Hep147_S196_R1_001.fastq.gz \
/nobackup1/kunyin71/Illumina_DNA_Reads/raw/Hep147_S196_R2_001.fastq.gz \
/nobackup1/kunyin71/Illumina_DNA_Reads/trimmed/Hep147_S196_R1_001_paired.fastq.gz \
/nobackup1/kunyin71/Illumina_DNA_Reads/trimmed/Hep147_S196_R1_001_unpaired.fastq.gz \
/nobackup1/kunyin71/Illumina_DNA_Reads/trimmed/Hep147_S196_R2_001_paired.fastq.gz \
/nobackup1/kunyin71/Illumina_DNA_Reads/trimmed/Hep147_S196_R2_001_unpaired.fastq.gz \
ILLUMINACLIP:/home/kunyin71/software/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36

```
### Quick GUI quality score
FastQC is a quick GUI quality score viewer. The quality score is defined as -10log<sub>10</sub>(P<sub>error</sub>).
 
```
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
fastqc read1.fastq.gz read2.fastq.gz ## it gives you a html file to view the quality score
```
## Sequence alignment
To begin with, you will need your filtered metagenomics data and the reference sequence. If you are downloading sequence from NCBI, please select Refseq option cause this gives you the original sequence. 
### Install Bowtie2
```
wget https://github.com/BenLangmead/bowtie2/archive/refs/heads/master.zip
make ## you will need to cd into the main folder
```
### Example script
This software compare all the short reads with the reference genome. In my case I need to get rid of all the sequences that match with the reference genome. So I use "--un-conc". It gave me two files with the shrimp DNA deleted. -S /dev/null here trash the default SAM file given by Bowtie2.

```
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 10-00:00:00
#SBATCH -p sched_mit_chisholm
#SBATCH -J skewer_Hep20_S195
#SBATCH --mem=250G

# Create a directory for the results
mkdir -p /home/kunyin71/reference_alignment/Hep20_S193

# Run Bowtie2 with the output directed to the new directory
/home/kunyin71/software/bowtie2/bowtie2 -x reference_index \
-1 /nobackup1/users/kunyin71/Illumina_DNA_Reads/trimmed/Hep20_S193-trimmed-pair$
-2 /nobackup1/users/kunyin71/Illumina_DNA_Reads/trimmed/Hep20_S193-trimmed-pair$
--un-conc /home/kunyin71/reference_alignment/Hep20_S193/unaligned_%.fastq \
-S /dev/null
```
### Bowtie2 output
Finally I get the genome fastq file and a log telling me how much of the reads turn out to be the shrimp DNA reads, as follows.

```
19190906 reads; of these:
  19190906 (100.00%) were paired; of these:
    11087731 (57.78%) aligned concordantly 0 times
    4633744 (24.15%) aligned concordantly exactly 1 time
    3469431 (18.08%) aligned concordantly >1 times
    ----
    11087731 pairs aligned concordantly 0 times; of these:
      1192334 (10.75%) aligned discordantly 1 time
    ----
    9895397 pairs aligned 0 times concordantly or discordantly; of these:
      19790794 mates make up the pairs; of these:
        14664140 (74.10%) aligned 0 times
        2624875 (13.26%) aligned exactly 1 time
        2501779 (12.64%) aligned >1 times
61.79% overall alignment rate
```

## Build contigs

### Install megahit
```
wget https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz

tar zvxf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
cd MEGAHIT-1.2.9-Linux-x86_64-static/bin/

./megahit --test  # run on a toy dataset
```

### Syntax, building contigs with filtered sequence
```
./megahit -1 MYREAD1.fq.gz -2 MYREAD2.fq.gz -o MYOUTPUTDIR
```
### Example script
The reason why we build contigs with the six samples processed together is that those samples comes from similar environement and have similar genes. Therefore, the software can identify some shared features among samples and improve the final accuracy of grouping genes into contigs.

```
#!/bin/bash
#SBATCH -N 1 ## 
#SBATCH -n 20 ## cores, chisholm have 10 nodes each with 250 GB ram and 20 cpus/cores
#SBATCH -t 10-00:00:00 ## time limit, max in Chisholm is 90 days
#SBATCH -p sched_mit_chisholm ## tell it which scheduler
#SBATCH -J skewer_SAMPLE ## job name, unique for each one
#SBATCH --mem=250G ## if mapping, go a little bit above the database size
#SBATCH --job-name=sixsample_metagenomics
#SBATCH --output=/home/kunyin71/build_contigs/sixsample_%j.out
#SBATCH --error=/home/kunyin71/build_contigs/sixsample_%j.err

/home/kunyin71/software/megahit/bin/megahit \
-1 /home/kunyin71/reference_alignment/Hep174_S194/unaligned_1.fastq,/home/kunyin71/reference_alignment/Hep147_S196/unaligned_1.fastq,/home/kunyin71/reference_alignment/Hep184_S197/unaligned_1.fastq,/home/kunyin71/reference_alignment/Hep225_S198/unaligned_1.fastq,/home/kunyin71/reference_alignment/Hep20_S193/unaligned_1.fastq,/home/kunyin71/reference_alignment/hep161/unaligned_1.fastq \
-2 /home/kunyin71/reference_alignment/Hep174_S194/unaligned_2.fastq,/home/kunyin71/reference_alignment/Hep147_S196/unaligned_2.fastq,/home/kunyin71/reference_alignment/Hep184_S197/unaligned_2.fastq,/home/kunyin71/reference_alignment/Hep225_S198/unaligned_2.fastq,/home/kunyin71/reference_alignment/Hep20_S193/unaligned_2.fastq,/home/kunyin71/reference_alignment/hep161/unaligned_2.fastq \
-o /home/kunyin71/build_contigs/results \
-t 20

```
### Gene prediction with prodigal

```
#!/bin/bash
#SBATCH -N 1 ## 
#SBATCH -n 20 ## cores, chisholm have 10 nodes each with 250 GB ram and 20 cpus/cores
#SBATCH -t 10-00:00:00 ## time limit, max in Chisholm is 90 days
#SBATCH -p sched_mit_chisholm ## tell it which scheduler
#SBATCH -J skewer_SAMPLE ## job name, unique for each one
#SBATCH --mem=250G ## if mapping, go a little bit above the database size
#SBATCH --job-name=single_%j
#SBATCH --output=single_%j.out
#SBATCH --error=single_%j.err

prodigal -i /nobackup1/kunyin71/contigs_single/Hep161_S195_results/final.contigs.fa -o hep161.genes -a hep161.faa -p meta

prodigal -i  -o hep174.genes -a hep161.faa -p meta
```

### Mapping reads with contigs

build reference library with the assembled contigs

```
#!/bin/bash
#SBATCH -N 1 ## 
#SBATCH -n 20 ## cores, chisholm have 10 nodes each with 250 GB ram and 20 cpus/cores
#SBATCH -t 10-00:00:00 ## time limit, max in Chisholm is 90 days
#SBATCH -p sched_mit_chisholm ## tell it which scheduler
#SBATCH -J skewer_SAMPLE ## job name, unique for each one
#SBATCH --mem=250G ## if mapping, go a little bit above the database size
#SBATCH --job-name=refer_back
#SBATCH --output=/nobackup1/kunyin71/map_contigs/refer_%j.out
#SBATCH --error=/nobackup1/kunyin71/map_contigs/refer_%j.err

/home/kunyin71/software/bowtie2/bowtie2-build /home/kunyin71/build_contigs/results/final.contigs.fa /nobackup1/kunyin71/map_contigs/reference_index

```
Irratative script

```
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 10-00:00:00
#SBATCH -p sched_mit_chisholm
#SBATCH -J bowtie2_alignment_${1}
#SBATCH --mem=250G
#SBATCH -o bowtie2_alignment_%j.out
#SBATCH -e bowtie2_alignment_%j.err

# Sample name from command line argument
SAMPLE_NAME=$1

# Define base directory for reads and results
READS_DIR=/nobackup1/kunyin71/Illumina_DNA_Reads/trimmed
REFERENCE_INDEX=/home/kunyin71/map_contigs/reference_index
RESULTS_DIR=/nobackup1/kunyin71/map_contigs

# Run Bowtie2
/home/kunyin71/software/bowtie2/bowtie2 -p 20 -x $REFERENCE_INDEX \
-1 $READS_DIR/${SAMPLE_NAME}-trimmed-pair1.fastq.gz \
-2 $READS_DIR/${SAMPLE_NAME}-trimmed-pair2.fastq.gz \
-S $RESULTS_DIR/${SAMPLE_NAME}_aligned.sam

```

call the above script

```
#!/bin/bash

# Array of sample names
SAMPLE_NAMES=("Hep174_S194" "Hep147_S196" "Hep184_S197" "Hep225_S198" "Hep20_S193" "Hep161_S195")

# Iterate over sample names and call the script for each
for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"
do
    sbatch run_bowtie2.sh $SAMPLE_NAME
done
```
### Convert .SAM file to .BAM file
Following the instructions, install samtools, bcftools and htslib from this website. https://www.htslib.org/download/



## Build metagenome-assembled genomes
Hereby are the top three most widely-used MAG assembling tools. The ideal case is that you have three set of results coming individually from the three tools and finally consolidate them together.
### Install CONCOCT

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create --prefix /home/kunyin71/software/conda_envs/concoct_env python=3 concoct
conda activate /home/kunyin71/software/conda_envs/concoct_env

```

### using CONCOCT
conda activate /home/kunyin71/software/conda_envs/concoct_env

```
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 10-00:00:00
#SBATCH -p sched_mit_chisholm
#SBATCH -J concoct
#SBATCH --mem=250G

SAMPLE_NAME=$1
DIR=/nobackup1/kunyin71/concoct/$SAMPLE_NAME
mkdir $DIR

/home/kunyin71/software/conda_envs/concoct_env/bin/cut_up_fasta.py /nobackup1/kunyin71/build_contigs/results/final.contigs.fa -c 10000 -o 0 --merge_last -b $DIR/contigs_10K.bed > $DIR/contigs_10K.fa

/home/kunyin71/software/conda_envs/concoct_env/bin/concoct_coverage_table.py $DIR/contigs_10K.bed /nobackup1/kunyin71/sam_to_bam/sorted_${SAMPLE_NAME}_aligned.bam > $DIR/coverage_table.tsv

concoct --composition_file $DIR/contigs_10K.fa --coverage_file $DIR/coverage_table.tsv -b $DIR/

/home/kunyin71/software/conda_envs/concoct_env/bin/merge_cutup_clustering.py $DIR/clustering_gt1000.csv > $DIR/clustering_merged.csv

mkdir $DIR/fasta_bins

/home/kunyin71/software/conda_envs/concoct_env/bin/extract_fasta_bins.py /nobackup1/kunyin71/build_contigs/results/final.contigs.fa $DIR/clustering_merged.csv --output_path $DIR/fasta_bins
```
### Install MaxBin
Download from https://sourceforge.net/projects/maxbin2/files/latest/download, cd into the src file, ```make```, ```cd ..```

Do ```./autobuild_auxiliary``` and you will be good. 

Call ```perl run_MaxBin.pl -h``` to see how to use the software.
### Install MetaBAT
You need to load a bunch of module from MIT public shared module list, including cmake, gcc, autoconf, and boost.

Specially, for boost, you may need to install it in your own directory. Detailed instructions are provided here. You might need to go through 5.1 Easy build and install. The tricky part here is the storage. Boost is a C++ library and contains 70,000+ scripts. Unfortunately, MIT engaging cluster is very vulnerable when being written a lot of small files in. It's going to take time.
```
https://www.boost.org/doc/libs/1_84_0/more/getting_started/unix-variants.html
```
You need to run the following commands in the main folder to have boost compiled to reach the minimum requirement of metabat.

```
./bootstrap.sh
./b2 --with-program_options --with-filesystem --with-system --with-graph --with-serialization --with-iostreams
```
For the rest of the steps, you should get the latest version from their website, but it should be something as follows.

Get the source code

```
wget https://bitbucket.org/berkeleylab/metabat/get/10d6f8779fca.zip

```

Compile it with something like

```
#!bash
#clean up old files
rm -f master.tar.gz
rm -f dev.tar.gz
rm -rf berkeleylab-metabat-*

#stable release version
wget https://bitbucket.org/berkeleylab/metabat/get/master.tar.gz
tar xzvf master.tar.gz
cd berkeleylab-metabat-*

# OR for the latest development version, uncomment these lines:
#wget https://bitbucket.org/berkeleylab/metabat/get/dev.tar.gz
#tar xzvf dev.tar.gz
#cd berkeleylab-metabat-*

#run the installation script
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/home/kunyin71/software/metabat -DBOOST_ROOT=home/kunyin71/software/boost/
make && make test && make install

```

I got a conflict between gcc and cmake. Need to compile either of them manually.

### MAG consolidation
Directly do
```
wget https://github.com/cmks/DAS_Tool/archive/refs/heads/master.zip
```

## Post-processing of metagenome-assembled genomes

### MAG evaluation
