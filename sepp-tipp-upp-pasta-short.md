``` bash
git clone git@github.com:smirarab/pasta.git 
git clone git@github.com:smirarab/sepp.git
git clone git@github.com:smirarab/sate-tools-mac.git
cd pasta/
sudo python setup.py develop 
pwd
run_pasta.py -h
cd ../sepp
python setup.py config -c
sudo python setup.py develop
run_sepp.py -h
which run_pasta.py
python setup.py upp -c
run_upp.py -h
wget http://www.cs.utexas.edu/~phylo/software/sepp/tipp.zip
curl -O http://www.cs.utexas.edu/~phylo/software/sepp/tipp.zip
unzip tipp.zip 
export REFERENCE=`pwd`/tipp
which blastn
export BLAST=`which blastn`
python setup.py tipp -c
run_tipp.py -c
run_tipp.py -h
cd ..
ls
mkdir seppRuns
cd seppRuns
ln -s ../sepp/test/unittest/data/mock/ .
run_sepp.py -t mock/rpsS/sate.tre -r mock/rpsS/sate.tre.RAxML_info -a mock/rpsS/sate.fasta -f mock/rpsS/rpsS.even.fas -o rpsS.out.default
pwd
cd ..
mkdir tippRuns
cd tippRuns/
ln -s ../sepp/test/ .
run_tipp.py -R pyrg -f test/unittest/data/mock/pyrg/pyrg.even.fas  -o output -P 30
head output_classification.txt 
mkdir profile
run_tipp_tool.py -g pyrg -a profile -o profile -p pyrg -i output_classification.txt -t 0.95
head profile/pyrg.classification 
head profile/abundance.species.csv 
cp test/unittest/data/mock/mixed/facs_simhc.short.fas .
run_abundance.py -f facs_simhc.short.fas -c ../sepp/.sepp/tipp.config -d out
head out/markers/abundance.species.csv
ls
head out/abundance.species.csv 
pwd
cd ../
mkdir uppRuns
cd uppRuns/
ln -s ../sepp/test/ . -s 
ln -s ../sepp/test/ .
run_upp.py -A 10 -B 1000 -M -1 -m amino -s test/unittest/data/upp_frag/amino.fas
pwd
ls -alt
head -n4 output_alignment_masked.fasta
pwd
cd ..
mkdir pastaRuns
cd pastaRuns/
ln -s ../pasta/data/ .
run_pasta.py -i data/small.fasta
cd ..
history > tutorial.md
```
