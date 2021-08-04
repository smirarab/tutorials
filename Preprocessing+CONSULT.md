
# CONSULT
<!-- Accurate contamination removal using locality-sensitive hashing-->

CONSULT is a tool for contamination removal from genomic sequencing reads. Relying on locality-sensitive hashing, CONSULT extracts *k*-mers from a query set and tests whether they fall within a user-specified hamming distance of *k*-mers in the reference dataset. It supports the inclusion of approximately 8 billion *k*-mers in its reference library, accommodating datasets with tens of thousands of microbial species.

The paper where we have described design of the algorithm and software architecture are now available on online (https://www.biorxiv.org/content/10.1101/2021.03.18.436035v1). <!-- (open access): -->
<!--  - [paper reference and doi][1] -->

* Summary data tables and scripts that we used during testing are available at https://github.com/noraracht/lsh_scripts.

* Raw data are deposited in https://github.com/noraracht/lsh_raw_data.

* Our custom CONSULT libraries constructed using different genomic reference sets:
    - [GTDB](https://tera-trees.com/data/consult/v1.0.0/all_nbrhood_kmers_k32_p3l2clmn7_K15-map2-171_gtdb.tar.gz)
    - [TOL](https://tera-trees.com/data/consult/v1.0.0/all_nbrhood_kmers_k32_p3l2clmn7_K15-map2-171_ToL.tar.gz)
    - [Bacterial/Archaeal Kraken](https://tera-trees.com/data/consult/v1.0.0/all_nbrhood_kmers_k32_p3l2clmn7_K15-map2-171_kraken.tar.gz)
    - [Mitochondrial CONSULT database](https://tera-trees.com/data/consult/v1.0.0/consult_mito_k32_p3l2clmn7_K15_tag2_v171.tar.gz)
    - [Plastid CONSULT database](https://tera-trees.com/data/consult/v1.0.0/consult_plastid_k32_p3l2clmn7_K15_tag2_v171.tar.gz)
  
  <!--  - [GTDB](https://drive.google.com/file/d/1MQJAXmZiTurumlZpvNoMLB0tKWGM_VE4/view?usp=sharing)-->
  <!--    - [TOL](https://drive.google.com/file/d/1sA9HFjWoU2jZ2vjd98pHVDEFRzOKMImk/view?usp=sharing)-->
  <!--    - [Bacterial/Archaeal Kraken](https://drive.google.com/file/d/1jeZB6b6aXl06BpPPsjM8oQA4xingJ1Dq/view?usp=sharing)-->
  <!--   - [Mitochondrial CONSULT database](https://drive.google.com/file/d/1mFD3dYFrJKqUkWlkRHbrQt-6eG-_K5vI/view?usp=sharing)-->

At the moment when using our libraries, library name needs to stay unchanged since library files are prefixed with library name. This will change in future releases to allow for more flexibility.
