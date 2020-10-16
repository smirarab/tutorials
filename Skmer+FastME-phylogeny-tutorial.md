# A tutorial on phylogenetic reconstruction using skimming data 


This tutorial walks you through the use of Skmer for computing distances (after data cleanup using Kraken and BBTools) between genome skims, FastMe for phylgoenetic reconstruction.

See [this tutorial](Skmer+APPLES+tutorial.md) for the closely related topic of phylogenetic placement. 


### Tools:

Tools used:

* [Skmer](https://github.com/shahab-sarmashghi/Skmer) for distance calculation from skims
* [FastME](http://www.atgc-montpellier.fr/fastme/binaries.php) for phylogeny reconstruction
* [BBTools](https://sourceforge.net/projects/bbmap) for read cleanup
* [Kraken-II](https://github.com/DerrickWood/kraken2.git) for filtering bacterial contamination

**System requirement:**

- Of the tools mentioned here, Kraken-II is practically only available on Linux servers with high enough memory and disk space. Others should work on MAC as well. If you cannot run the Kraken step, you can skip it. 


**NOTE:** 
Of all of our steps, Kraken-II filtering is the only "trick" one. It requires powerful machines, Linux, and lots of space. So you may want to skip it in your first try. 

# Installations

* Instructions shown below work on Linux. Some changes are needed for other systems. See comments.

### Install tools

~~~bash
mkdir tutorial
cd tutorial 

conda create --name tutorial
conda activate tutorial

conda config --add channels bioconda

### Instal Skmer
conda install skmer
skmer -h


### Install FastME (to get backbone trees)
wget http://www.atgc-montpellier.fr/download/sources/fastme/fastme-2.1.5.tar.gz
tar xvfz fastme-2.1.5.tar.gz
chmod +x fastme-2.1.5/binaries/fastme-2.1.5-linux64 ## Change "linux64" at the end if using other platforms (osx or windows).
./fastme-2.1.5/binaries/fastme-2.1.5-linux64 -h


### download bbtools and extract archive
wget -O bbmap.tar.gz  https://sourceforge.net/projects/bbmap/files/BBMap_38.87.tar.gz/download
tar -zxvf bbmap.tar.gz
rm bbmap.tar.gz

### download bbtools preprocessing pipeline
wget https://raw.githubusercontent.com/balabanmetin/tools/master/bbmap_pipeline.sh
chmod +x ./bbmap_pipeline.sh

### clone kraken2 source
git clone https://github.com/DerrickWood/kraken2.git

### install kraken2 software
cd kraken2 && ./install_kraken2.sh . && cd ..


### download a small script that converts skmer output format to FastME input format
wget https://raw.githubusercontent.com/balabanmetin/LSEdiag/master/tsv_to_phymat.sh
~~~




### Obtain large datasets





~~~bash
### Obtain yeast genomes as test case
wget https://github.com/balabanmetin/yeast-genomes/raw/master/yeast-genomes.tar.bz2
tar xvfj yeast-genomes.tar.bz2
rm yeast-genomes.tar.bz2

du -sm genomes/*/*

head genomes/Saccharomyces_kudriavzevii/GCA_900682665.1_SKCA111_genomic.fna


### Obtain genome skims simulated based on these genomes
wget https://github.com/balabanmetin/yeast-genomes/raw/master/yeast-genome-skims.tar.bz2
tar -jxvf yeast-genome-skims.tar.bz2
rm yeast-genome-skims.tar.bz2


### We want to now download a bunch of files from this google drive link: 
###     https://drive.google.com/drive/folders/1vDfSY2o9-zXPY5a0E64oWjqRaqPyCrsL?usp=sharing
# I like to do do things in commandline, so I need to download a tool for it.
pip install gdown
mkdir kraken_db_standard_k35l31s7
cd kraken_db_standard_k35l31s7
gdown  "https://drive.google.com/uc?id=1389YOsQ9skqecyRImry_vogl9TVqcfBc"
gdown  "https://drive.google.com/uc?id=1ved2D016ZxmpoUNm3mDIp8O4FCs05f_G"
# This is a big database (40 GB) and it may take a while to finish.
gdown  "https://drive.google.com/uc?id=1bzgxAHrOjSchF3B6nOaaQnYMP19VOLlG"
cd ..
~~~



# Analyses using genome skims

Before you start, note that we have a file that gives the list of species names. We will use this file to make sure we only work on genomes we want to analyze. Note that are `.fq` read files are named consistently with this list. In your datasets, if you have different naming conventions, please make a file that has the names appearing in `.fq` files. 

~~~bash
cat genomes/nonhybrids.txt
~~~

## Cleanup skims using BBtools
~~~bash

# This script dedupblicates, removes primers, and merges paired-end reads
for x in `cat genomes/nonhybrids.txt`; do 
        ./bbmap_pipeline.sh reads/${x}1.fq reads/${x}2.fq preprocessed/${x}.fq
done # We only choose non-hybrid species
~~~

## Kraken filtering

* If you have Kraken-II, run:
* 
~~~bash
mkdir -p {classified,unclassified,reports}
~~~

* If you don't have Kraken-II, let's assume for now we have no contamination and run:
* 
~~~bash
mkdir -p {unclassified}
for x in `cat genomes/nonhybrids.txt`; do 
        cp preprocessed/${x}.fq unclassified/
done 
~~~

## Compute distances  using genome skims

~~~bash
# unclassified directory contains fastq files without contamination.
# Now, run Skmer on this input directory.
skmer reference -t unclassified/
# This command will sketch each input file and will also compare all pairs to compute their distance
# On my machine it takes 4-5 minutes
# `-t` Tells is to transform distances using JC69 transformation

# Look at distances computed
cat ref-dist-mat.txt

# Also look at genome length and estimated coverage
head library/*/*.dat
~~~

## Infer backbone tree  using FastME
~~~bash
# Reformat distance file
bash tsv_to_phymat.sh ref-dist-mat.txt ref-dist-mat.phy
cat ref-dist-mat.phy


fastme-2.1.5/binaries/fastme-2.1.5-linux64  -i ref-dist-mat.phy  -o skim-phylogeny.tre
~~~

# Alignment-free analyses using full genomes

Skmer overwrites existing files with a fixed output name. Thus, to avoid loosing our library based on skims, we rename the outputs (sorry, this should be cleaner in future)/ 

~~~bash
mv library library-skims
mv ref-dist-mat.txt ref-dist-mat-skims.txt
~~~


## Run Skmer on example genome data 

~~~bash
# First, prepare a directory with all the input files
mkdir nonhybrids
for x in `cat genomes/nonhybrids.txt`; do  
	cp genomes/$x/*fna  nonhybrids/$x.fna; 
done # We only choose non-hybrid species

# Now, run Skmer on this input directory.
skmer reference -t nonhybrids/
# This command will sketch each input file and will also compare all pairs to compute their distance
# On my machine it takes 4-5 minutes
# `-t` Tells is to transform distances using JC69 transformation
~~~

## Infer backbone tree  using FastME
~~~bash
bash tsv_to_phymat.sh ref-dist-mat.txt ref-dist-mat.phy

# Infer backbone from scratch
fastme-2.1.5/binaries/fastme-2.1.5-osx  -i ref-dist-mat.phy  -o genome-phylogeny.tre
~~~


If you have a backbone tree, compute branch lengths on it, you can use FastME as well.. 

~~~bash
echo "((((Saccharomyces_paradoxus,(Saccharomyces_jurei,Saccharomyces_mikatae)),Saccharomyces_kudriavzevii),Saccharomyces_arboricola),(Saccharomyces_eubayanus,Saccharomyces_uvarum));" > backbone.tre
# Recompute branch lengths:
fastme-2.1.5/binaries/fastme-2.1.5-osx  -i ref-dist-mat.phy  -o backbone-fastme-2.tre -u backbone-fastme.tre
~~~




# Papers:


* S. Sarmashghi, K. Bohmann, M. T. P Gilbert, V. Bafna, and S. Mirarab. “Skmer: Assembly-Free and Alignment-Free Sample Identification Using Genome Skims.” Genome Biology Vol. 20, no. 1 (2019): pp. 34. doi:10.1186/s13059-019-1632-4.

Other relevant papers:

* K. Bohmann, S. Mirarab, V. Bafna, and M. T. P. Gilbert. “Beyond DNA Barcoding: The Unrealized Potential of Genome Skim Data in Sample Identification.” Molecular Ecology, (2020), pp. mec.15507. doi:10.1111/mec.15507.
* E. Rachtman, M. Balaban, V. Bafna, and S. Mirarab. “The Impact of Contaminants on the Accuracy of Genome Skimming and the Effectiveness of Exclusion Read Filters.” Molecular Ecology Resources Vol. 20, no. 3 (2020): pp. 1755-0998.13135. doi:10.1111/1755-0998.13135.
* M. Balaban, S. Sarmashghi, and S. Mirarab. “APPLES: Scalable Distance-Based Phylogenetic Placement with or without Alignments.” Edited by David Posada. Systematic Biology Vol. 69, no. 3 (2020): pp. 566–78. doi:10.1093/sysbio/syz063.
* M. Balaban, and S. Mirarab. “Phylogenetic Double Placement of Mixed Samples.” Bioinformatics Vol. 36, no. Supplement_1 (2020): pp. i335–43. doi:10.1093/bioinformatics/btaa489.

