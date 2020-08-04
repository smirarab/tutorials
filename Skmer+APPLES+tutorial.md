# Skmer+APPLES+MISA tutorial

This toturial walks you through the use of Skmer for computing distances between genome skims, APPLES for phylgoenetic placement, and MISA for phylgoenetic placement of mixed samples (of two species). 

### Tools:

Main tools:

* [Skmer](https://github.com/shahab-sarmashghi/Skmer)
* [APPLES](https://github.com/balabanmetin/apples)
* [MISA](https://github.com/balabanmetin/misa)


Other tools we will use:

* [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/)
* [FastME](http://www.atgc-montpellier.fr/fastme/binaries.php)
* [FastTree-II](http://www.microbesonline.org/fasttree)
* [guppy](https://matsen.github.io/pplacer/generated_rst/guppy.html)

### Papers:


* S. Sarmashghi, K. Bohmann, M. T. P Gilbert, V. Bafna, and S. Mirarab. “Skmer: Assembly-Free and Alignment-Free Sample Identification Using Genome Skims.” Genome Biology Vol. 20, no. 1 (2019): pp. 34. doi:10.1186/s13059-019-1632-4.
* M. Balaban, S. Sarmashghi, and S. Mirarab. “APPLES: Scalable Distance-Based Phylogenetic Placement with or without Alignments.” Edited by David Posada. Systematic Biology Vol. 69, no. 3 (2020): pp. 566–78. doi:10.1093/sysbio/syz063.
* M. Balaban, and S. Mirarab. “Phylogenetic Double Placement of Mixed Samples.” Bioinformatics Vol. 36, no. Supplement_1 (2020): pp. i335–43. doi:10.1093/bioinformatics/btaa489.
* K. Bohmann, S. Mirarab, V. Bafna, and M. T. P. Gilbert. “Beyond DNA Barcoding: The Unrealized Potential of Genome Skim Data in Sample Identification.” Molecular Ecology, (2020), pp. mec.15507. doi:10.1111/mec.15507.
* E. Rachtman, M. Balaban, V. Bafna, and S. Mirarab. “The Impact of Contaminants on the Accuracy of Genome Skimming and the Effectiveness of Exclusion Read Filters.” Molecular Ecology Resources Vol. 20, no. 3 (2020): pp. 1755-0998.13135. doi:10.1111/1755-0998.13135.

# Installations

Instructions shown below work on MAC OS Catalina. Some changes are needed for Linux and Windows. See comments.

### Install main tools

~~~bash
conda create --name tutorial
conda activate tutorial

conda config --add channels bioconda

### Instal Skmer
conda install skmer
skmer -h

### Instal APPLES
python -m pip install -U apples
run_apples.py -h
python -m pip  list |grep apples
### If you have versions older than 1.3.0, you may need to updating using:
python -m pip install --upgrade apples

### Install MISA
python -m pip install misa
run_misa.py -h
~~~

### Install other methods
~~~bash
### Install FastME (to get backbone trees)
wget http://www.atgc-montpellier.fr/download/sources/fastme/fastme-2.1.5.tar.gz
tar xvfz fastme-2.1.5.tar.gz
chmod +x fastme-2.1.5/binaries/fastme-2.1.5-osx ## Change "osx" at the end if using other platforms (inux or windows).
./fastme-2.1.5/binaries/fastme-2.1.5-osx -h

### install fasttree
wget http://www.microbesonline.org/fasttree/FastTree.c
gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm
### Note: for linux, you can find binary files on the FastTree website
./FastTree -h

### Guppy for working with .jplace files
wget https://github.com/smirarab/sepp/raw/master/tools/bundled/Darwin/guppy
# wget https://github.com/smirarab/sepp/raw/master/tools/bundled/Linux/guppy-64 for linux
chmod +x guppy
./guppy --version

### (Optional) If you want to simulate skimming, Install ART
# The link below is for MAC. For other platforms, see https://www.niehs.nih.gov/research/resources/software/biostatistics/art/
    # For Linux, it's https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05macos64.tgz
tar xvfj artbinmountrainier2016.06.05macos64.tgz
art_bin_MountRainier/art_illumina -h
~~~

### Obtain large datasets

~~~bash
### Obtain yeast genomes as test case
wget https://github.com/balabanmetin/yeast-genomes/raw/master/yeast-genomes.tar.bz2
tar xvfj yeast-genomes.tar.bz2

du -sm genomes/*/*

head genomes/Saccharomyces_kudriavzevii/GCA_900682665.1_SKCA111_genomic.fna
~~~

# Alignment-free analyses: full genomes

## Run Skmer on example genome data 
~~~bash
# First, prepare a directory with all the input files
mkdir nonhybrids
for x in `cat genomes/nonhybrids.txt`; do  
	cp genomes/$x/*fna  nonhybrids/$x.fna; 
done # We only choose non-hybrid species
# Then, set aside the genome S. cerevisiae as query
mkdir nonhybrids-query 
mv nonhybrids/Saccharomyces_cerevisiae.fna nonhybrids-query/


# Now, run Skmer on this input directory.
skmer reference -t nonhybrids/
# This command will sketch each input file and will also compare all pairs to compute their distance
# On my machine it takes 4-5 minutes
# `-t` Tells is to transform distances using JC69 transformation
~~~

## Infer backbone tree with branch lengths using FastME
~~~bash
### download a small script that converts skmer output format to FastME input format
wget https://raw.githubusercontent.com/balabanmetin/LSEdiag/master/tsv_to_phymat.sh
bash tsv_to_phymat.sh ref-dist-mat.txt ref-dist-mat.phy

# Infer backbone from scratch
fastme-2.1.5/binaries/fastme-2.1.5-osx  -i ref-dist-mat.phy  -o backbone-fastme.tre

# If you have a backbone tree, compute branch lengths on it
echo "((((Saccharomyces_paradoxus,(Saccharomyces_jurei,Saccharomyces_mikatae)),Saccharomyces_kudriavzevii),Saccharomyces_arboricola),(Saccharomyces_eubayanus,Saccharomyces_uvarum));" > backbone.tre
# Recompute branch lengths:
fastme-2.1.5/binaries/fastme-2.1.5-osx  -i ref-dist-mat.phy  -o backbone-fastme-2.tre -u backbone-fastme.tre

~~~


## Run actual placement
~~~bash
# A small script again to convert formats
wget https://raw.githubusercontent.com/balabanmetin/misa/master/scripts/convert_to_tsv.sh

# Ask skmer to compute distances from query to references
skmer query -t nonhybrids-query/Saccharomyces_cerevisiae.fna library/ # Use -a if you want the query to be added
# Convert to right format using script above
bash convert_to_tsv.sh dist-saccharomyces_cerevisiae.txt > dist-saccharomyces_cerevisiae.tsv
run_apples.py -t backbone-fastme.tre -d dist-saccharomyces_cerevisiae.tsv -o cerevisiae.jplace

# Turn output to newick
./guppy tog cerevisiae.jplace
~~~


# Alignment-free analyses: genome skims

## Run ART to simulate genome skims
~~~bash
# First, use ART to create genome skims from your genomes at 2X coverage
mkdir skims skims/nonhybrids skims/nonhybrids-query/
ls nonhybrids*/*fna|xargs -n1 -I@ art_bin_MountRainier/art_illumina  -i @  -l 100 -f 2 -o skims/@
cd skims
~~~

## Repeat analyses using genome skims
~~~bash
# Build library of skims
skmer reference -t nonhybrids/

# Compute distance of query to backbone sequence
skmer query -t nonhybrids-query/Saccharomyces_cerevisiae.fna.fq library/ 
sed -i -e "s/.fna//g" dist-saccharomyces_cerevisiae.fna.txt # Remove extra .fna from species names.
# Convert file format
bash ../convert_to_tsv.sh dist-saccharomyces_cerevisiae.fna.txt > dist-saccharomyces_cerevisiae.tsv


# Run APPLES
run_apples.py -t ../backbone-fastme.tre -d dist-saccharomyces_cerevisiae.tsv -o cerevisiae.jplace

# Turn output to newick
./guppy tog cerevisiae.jplace

# Compare to the genome-based run

cd ..
~~~



# Using MISA for mixed genome skim analyses 
~~~bash
mkdir mix-query
cp genomes/Saccharomyces_pastorianus/GCA_001515485.2_Saccharomyces_pastorianus_Weihenstephan_34_70_chromosomes_assembly_1.0_genomic.fna mix-query/Saccharomyces_pastorianus.fna

# Run Skmer
skmer query -t mix-query/Saccharomyces_pastorianus.fna library/
# Convert output to .tsv file
bash convert_to_tsv.sh dist-saccharomyces_pastorianus.txt > dist-saccharomyces_pastorianus.tsv

# Run MISA for phylogenetic double placemet
run_misa.py -d dist-saccharomyces_pastorianus.tsv -t backbone-fastme.tre -o mixed-output.jplace


# Check the output versus correct mixture:
./guppy tog mixed-output.jplace
cat genomes/Saccharomyces_pastorianus/things.txt
~~~

How about contamination?
 
* See Rachtman et al, 2019. Suggest using Kraken as one possible solution.
* Stay tuned ...

# Using APPLES on aligned sequences

~~~bash
# Get example input data
wget https://raw.githubusercontent.com/balabanmetin/apples/master/data/ref.fa
wget https://raw.githubusercontent.com/balabanmetin/apples/master/data/query.fa
wget https://raw.githubusercontent.com/balabanmetin/apples/master/data/backbone.nwk
sed -i -e "s/>L/>QU/g"  query.fa # Change query names to become easy to distinguish

# Recompute branch lengths on the backbone tree 
./FastTree  -nosupport -nt -nome -noml -log tree.log -intree backbone.nwk < ref.fa > minimum_evo_backbone.nwk

run_apples.py -s ref.fa -q query.fa -t minimum_evo_backbone.nwk -o aligned.jplace

# Turn output to newick
./guppy tog aligned.jplace
~~~
