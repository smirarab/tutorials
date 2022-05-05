Adopted from:

* [this tutorial](https://github.com/smirarab/tutorials/blob/master/Skmer-APPLES-tutorial.md)
* [and this tutorial](https://github.com/smirarab/tutorials/blob/master/Skmer%2BFastME-phylogeny-tutorial.md)

## Setup:

#### 1. Install Skmer
Requirements:

1. Conda with python3

~~~bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install skmer

skmer -h

~~~

### 2. Install APPLES-2

Requirements:

1. Pip with python3 (sometimes called pip3)

~~~bash
python -m pip install -U apples
run_apples.py -h
python -m pip  list |grep apples
### If you have versions older than 1.3.0, you may need to updating using:
python -m pip install --upgrade apples
~~~

#### 3. Install helper scripts and other tools

Instructions shown below work on **MAC OS** Catalina. Some changes are needed for Linux and Windows. See comments.

~~~bash
### Install fastme
wget http://www.atgc-montpellier.fr/download/sources/fastme/fastme-2.1.5.tar.gz
tar xvfz fastme-2.1.5.tar.gz
chmod +x fastme-2.1.5/binaries/fastme-2.1.5-linux64 ## Change to "osx" at the end if using MAC.
./fastme-2.1.5/binaries/fastme-2.1.5-linux64 -h


### Scripts to change data file formats:
wget https://raw.githubusercontent.com/balabanmetin/LSEdiag/master/tsv_to_phymat.sh
wget https://raw.githubusercontent.com/balabanmetin/misa/master/scripts/convert_to_tsv.sh
chmod +x tsv_to_phymat.sh convert_to_tsv.sh

### Guppy for working with .jplace files (not essential)
wget https://github.com/smirarab/sepp/raw/master/tools/bundled/Linux/guppy-64
# wget https://github.com/smirarab/sepp/raw/master/tools/bundled/Darwin/guppy for MacOS
chmod +x guppy
./guppy --version
~~~

These are useful but not needed right now. 

~~~bash
### install fasttree
wget http://www.microbesonline.org/fasttree/FastTree.c
gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm
### Note: for linux, you can find binary files on the FastTree website
./FastTree -h


### (Optional) If you want to simulate skimming, Install ART
# The link below is for MAC. For other platforms, see https://www.niehs.nih.gov/research/resources/software/biostatistics/art/
    # For Linux, it's https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05macos64.tgz
tar xvfj artbinmountrainier2016.06.05macos64.tgz
art_bin_MountRainier/art_illumina -h
~~~

#### 3. Get Data

Downloaded data from [here](https://drive.google.com/drive/folders/1r9IaMt392AY9UYSQwVSoTEhZB63Ar-Qz) onto **a folder called `skmer`**.

## Analyses of data

### Make a reference library

1. First, let's build the reference library.

    ~~~bash
	  skmer reference -t skims/
    ~~~
    * For large input files like these, you need ~5GB of memory per CPU core you use. You can use `-p` to adjust the number of cores
    * This can take a good 20 minutes on my machine.
2. Then, let's inspect it:

    ~~~bash
    # This is where the library is saved
	  ls library 
	  
	 # We can also look at estimates of coverage and error rate.
	 # These are kept in .dat files
	 
	 cat library/*/*.dat
	 
	 # Or grep all of them
	 grep coverage library/*/*dat

	 grep genome_length library/*/*dat
    ~~~
3. Note that sometimes Skmer fails to estimate coverage and error. In this dataset, several files are odd in that they have only 40 reads. These are to be discarded. To discard, we simply remove them. 

    ~~~bash
	 du -sk skims/*
	 # Note how some files are super small. 
	 
	 mkdir failedruns
	 # To remove, simply move the directories out of library
	 mv library/FSFRE3B597_1_S119_R1_001 library/FSFRE3B597_1_S119_R2_001 library/FSVIO6B764_1_S117_R1_001 library/FSVIO6B764_1_S117_R2_001 library/NLLARSB678_S114_R1_001 library/NLLARSB678_S114_R2_001 failedruns/
    ~~~
4. There is already a distance matrix that you can investigate:

   ~~~bash
     column -t ref-dist-mat.txt
   ~~~ 
	Note, however, that this file now includes those poor samples. 
5. Let's recreate the distance matrix with the new library. 

   ~~~bash
     skmer distance -t library 
     column -t ref-dist-mat.txt
   ~~~ 
   * The `-t` option simply instructs it add the Jukes-Cantor (69) phylogenetic correction
6. Look at distances in R

   ~~~bash
    echo "m=as.matrix(read.csv('ref-dist-mat.txt',sep='\t',row.names=1)); pdf('fig.pdf',width=12,height=9);plot(hclust(as.dist(m))); heatmap(m, scale = 'none');dev.off()"|R --vanilla
   ~~~ 
   Download and loot at the file called `fig.pdf`.

### Build a tree

1. Infer the tree

    ~~~bash
    # Convert data format
    bash tsv_to_phymat.sh ref-dist-mat.txt ref-dist-mat.phy
    
    # Run FastME (replace -linux with -osx for MAC)
    fastme-2.1.5/binaries/fastme-2.1.5-linux64  -i ref-dist-mat.phy  -o backbone.tre
    ~~~
    
2. Let's look at the tree:

   ~~~bash
   cat backbone.tre
   
   # If you have newick utilities:
   nw_reroot backbone.tre |nw_display -
			    | CEJES6B739 S115 R1 001
			    |
			    | CCALV3B613 S116 R2 001
			    |
			    |                               | |it 1938 S10 R2 001
			    |     +---------------------------+
			    |-----+                           +-+ Dit 1938 S10 R1 001
	 +------------------+     |
	 |                  | CCAL|3B613 S116 R1 001
	 |                  |
	=|                 ||CEJES6B739 S115 R2 001
	 |
	 |                 ||FCBDE6B758 S118 R1 001
	 +------------------+
			    | FCBDE6B758 S118 R2 001

	 |----------------|---------------|----------------|-----
	 0            0.001           0.002            0.003
	 substitutions/site
   ~~~
  
  ### Add a query onto the tree
  
  1. Download some new sample; we grabbed from [here](https://drive.google.com/drive/folders/1j4aRRs5-q1ZoQUg4a6NMvDw09Ux7F2RK).
     ~~~bash
     # compute distances from query to library
     # -a tells the skmer to add the query to the library
     skmer query -t -a ditGra1_S163_R1_001.fastq library
     
     # look at those distances
     column -t dist-ditgra1_s163_r1_001.txt
     
     ditGra1_S163_R1_001
     Dit_1938_S10_R1_001     0.004247606221784569
     Dit_1938_S10_R2_001     0.004262008940587567
     CCALV3B613_S116_R2_001  0.005809290678303558
     CCALV3B613_S116_R1_001  0.005967056863191565
     FCBDE6B758_S118_R2_001  0.0062520119487337895
     CEJES6B739_S115_R2_001  0.006313259420129632
     FCBDE6B758_S118_R1_001  0.006453031543366838
     CEJES6B739_S115_R1_001  0.006609567108046734
     ~~~
  
  2. Use APPLES to add query to existing tree

     ~~~bash
     # Another format conversion
     bash convert_to_tsv.sh dist-ditgra1_s163_r1_001.txt > dist-ditgra1_s163_r1_001.tsv
     
     # Now, run APPLES
     run_apples.py -t backbone.tre -d dist-ditgra1_s163_r1_001.tsv -o s163.jplace
     ~~~
     
  3. Inspect the output
     ~~~bash
     # First, look at the file
     cat s163.jplace
     
     # best to covert to newick
     ./guppy-64 tog s163.jplace
     
     # and display newick tree
     cat s163.tog.tre|nw_reroot -| nw_display -
		 +-----------------+ ditGra1 S163 R1 001
		 |
		 |                   ||Dit 1938 S10 R2 001
		=|                 +--+
		 |                 |  ++ Dit 1938 S10 R1 001
		 |                 |
		 +-----------------+         | C|ALV3B613 S116 R1 001
						   |            |
						   |            |  | CEJES6B739 S115 R1 001
						   +------------+  |
										|  | CCALV3B613 S116 R2 001
										+--+
										   | CEJES6B739 S115 R2 001
										   |
										   |                    | FCBDE6B758 S118 R1 001
										   +--------------------+
																| FCBDE6B758 S118 R2 001
		 |--------|---------|--------|--------|--------|---------
		 0    0.001     0.002    0.003    0.004    0.005
		 substitutions/site
     ~~~ 
     
  4. You could also redo the tree

	  ~~~bash
	  # Recompute distance matrix
	  skmer distance -t library
	  
	  # Look at distances
	  column -t ref-dist-mat.txt
	  
	  # Redo the figure
	  echo "m=as.matrix(read.csv('ref-dist-mat.txt',sep='\t',row.names=1)); pdf('fig2.pdf',width=12,height=9);plot(hclust(as.dist(m))); heatmap(m, scale = 'none');dev.off()"|R --vanilla
	  
	  #   Download and loot at the file called `fig2.pdf`.
	  
	  # Convert data format
	  bash tsv_to_phymat.sh ref-dist-mat.txt ref-dist-mat.phy
	    
	  # Run FastMEredo to redo the tree 
	  fastme-2.1.5/binaries/fastme-2.1.5-linux64  -i ref-dist-mat.phy  -o backbone.tre
	
	  # And display it:
	  nw_reroot backbone.tre |nw_display -
		 +-----------------+ ditGra1 S163 R1 001
		 |
		 |                                  | CEJES6B739 S115 R1 001
		 |                                  |
		 |                                  | CCALV3B613 S116 R2 001
		 |                                  |
		=|                               +--+                   | FCBDE6B758 S118 R2 001
		 |                               |  |-------------------+
		 |                               |  |                   | FCBDE6B758 S118 R1 001
		 |                 +-------------+  |
		 |                 |             |  | CEJES6B739 S115 R2 001
		 |                 |             |
		 +-----------------+          | C|ALV3B613 S116 R1 001
						   |
						   |||Dit 1938 S10 R2 001
						   +-+
							 ++ Dit 1938 S10 R1 001

		 |-----------------|-----------------|-----------------|-
		 0             0.002             0.004             0.006
		 substitutions/site
	  ~~~
  
 ### Add another query onto the tree 
 
 Repeat these steps with `ditGra2_S164_R2_001.fastq`. 
 
 ~~~bash
 skmer query -t -a ditGra2_S164_R2_001.fastq library
 
 skmer distance -t library

 # Look at distances
 column -t ref-dist-mat.txt
 
 # Update figure
 echo "m=as.matrix(read.csv('ref-dist-mat.txt',sep='\t',row.names=1)); pdf('fig3.pdf',width=12,height=9);plot(hclust(as.dist(m))); heatmap(m, scale = 'none');dev.off()"|R --vanilla

 # Download and inspect fig3.pdf
 
 bash tsv_to_phymat.sh ref-dist-mat.txt ref-dist-mat.phy

 # Run FastME (replace -linux with -osx for MAC)
  fastme-2.1.5/binaries/fastme-2.1.5-linux64  -i ref-dist-mat.phy  -o backbone.tre
 
 # nw_reroot backbone.tre |nw_display -
                                  | CEJES6B739 S115 R1 001
                                  |
                                  | CCALV3B613 S116 R2 001
                                  |
                                  | CEJES6B739 S115 R2 001
                                  |
                       +----------+                  | |it 1938 S10 R2 001
                       |          |  +-----------------+
                       |          +--+                 ++ Dit 1938 S10 R1 001
 +---------------------+             |
 |                     |          | C|ALV3B613 S116 R1 001
 |                     |
 |                     |               | FCBDE6B758 S118 R2 001
=|                     +---------------+
 |                                     | FCBDE6B758 S118 R1 001
 |
 |                     +------+ ditGra1 S163 R1 001
 +---------------------+
                    | d|tGra2 S164 R2 001

 |----------|-----------|----------|----------|----------
 0      0.001       0.002      0.003      0.004
 substitutions/site
 ~~~
  
 ## Should we have merged? 
 
 1. Install our merging scripts and tools:
 
     ~~~bash
	 ### download bbtools and extract archive
	 wget -O bbmap.tar.gz  https://sourceforge.net/projects/bbmap/files/BBMap_38.87.tar.gz/download
	 tar -zxvf bbmap.tar.gz
	 rm bbmap.tar.gz
	
	 ### download bbtools preprocessing pipeline
	 wget https://raw.githubusercontent.com/balabanmetin/tools/master/bbmap_pipeline.sh
	 chmod +x ./bbmap_pipeline.sh
    ~~~
    
 2. Add the OmniC samples to everything else:

     ~~~bash
     mv  ditGra2_S164* skims/
     ~~~
   
 2. Remove adaptors and Merge reads:

	 ~~~bash
	 for x in skims/*R1*; do ./bbmap_pipeline.sh $x ${x/R1/R2} ${x/R1/both}; done
	 
	 mv skims/*both* merged/
	 
	  mv library library-unmerged
	 ~~~
 
 3. Repeat analyses from the start, 
 
	 ~~~bash
	 skmer reference merged/
	 # This time it takes half the time roughly 
	 
	 column -t ref-dist-mat.txt 
	 
	 # move the problematic ones
	 mv library/FSFRE3B597_1_S119_both_001 library/FSVIO6B764_1_S117_both_001/ library/NLLARSB678_S114_both_001/ failedruns/
	 
	 skmer distance -t library
	  
	 # Look at distances
     column -t ref-dist-mat.txt
	 sample                    CCALV3B613_S116_both_001  CEJES6B739_S115_both_001  Dit_1938_S10_both_001  FCBDE6B758_S118_both_001  ditGra1_S163_both_001  ditGra2_S164_both_001
	 CCALV3B613_S116_both_001  0.0                       0                         0.0006931091362385231  0.0021092287912777923     0.003094254134256472   0.0027175476124217737
	 CEJES6B739_S115_both_001  0                         0.0                       0.0010941671840699913  0.002256308751831658      0.0035038030984755227  0.002948968474035177
	 Dit_1938_S10_both_001     0.0006931091362385231     0.0010941671840699913     0.0                    0.003814417288488263      0.003019312482994785   0.0036312702142234955
	 FCBDE6B758_S118_both_001  0.0021092287912777923     0.002256308751831658      0.003814417288488263   0.0                       0.00305936348863708    0.002721788925570516
	 ditGra1_S163_both_001     0.003094254134256472      0.0035038030984755227     0.003019312482994785   0.00305936348863708       0.0                    0
	 ditGra2_S164_both_001     0.0027175476124217737     0.002948968474035177      0.0036312702142234955  0.002721788925570516      0                      0.0
     
	 # Update figure
     echo "m=as.matrix(read.csv('ref-dist-mat.txt',sep='\t',row.names=1)); pdf('fig4.pdf',width=12,height=9);plot(hclust(as.dist(m))); heatmap(m, scale = 'none');dev.off()"|R --vanilla

     # Download and inspect fig4.pdf
 
     bash tsv_to_phymat.sh ref-dist-mat.txt ref-dist-mat.phy

     # Run FastME (replace -linux with -osx for MAC)
      fastme-2.1.5/binaries/fastme-2.1.5-linux64  -i ref-dist-mat.phy  -o backbone-merged.tre
 
     nw_reroot backbone-merged.tre |nw_display -
											+-+ CEJES6B739 S115 both 001
											|
						+-------------------+ +---------------+ Dit 1938 S10 both 001
						|                   +-+
		 +--------------+                 | CC|LV3B613 S116 both 001
		 |              |
		=|              +---------------------+ FCBDE6B758 S118 both 001
		 |
		 |              ++ ditGra1 S163 both 001
		 +--------------+
					  | |itGra2 S164 both 001

		 |--------|--------|--------|-------|--------|---------
		 0   0.0005    0.001   0.0015   0.002   0.0025
		 substitutions/site
	 ~~~
  
      **Note how OmniC results look much better now.**
