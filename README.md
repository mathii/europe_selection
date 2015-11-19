#Eight thousand years of natural selection in Europe

This is the code used for the selection analyis from the paper "Eight thousand years of natural selection in Europe" (Mathieson et al. 2015 - http://dx.doi.org/10.1101/016477). 

The commands in analysis/ANALYSIS should replicate the analysis and allow you to generate the figures shown in the paper. This includes calls to submit batch jobs using bsub, which you should replace with the appropriate calls for your own system, or just run jobs locally, which may be slow. 

Unfortunately the paths are hard coded, so it will be easier if you name this directory ~/selection/code and put the data in the directory ~/data/v8/  with two subdirectories; use/v81kg_europe2names.{snp,ind,geno} for the ancient genotype data merged with 1000 genomes; and reads/jj2.chrXX.readcounts.txt.gz for the read count data. The read count data has 4 columns: rsID, SampleID, REF, ALT. So if there are N samples and M SNPs, then there are N rows per SNP: NM rows intotal. 

The ancient genotype data is available from the [Reich lab website](http://genetics.med.harvard.edu/reich/Reich_Lab/Datasets.html). To get the read counts, aligned BAMs are available from the [European Nucleotide Archive](http://www.ebi.ac.uk/ena) under accession PRJEB11450.

The analysis will create various subdirectories in the ~/selection/ directory. If you have to put the data somewhere else, you'll have to change all the paths in the code manually. 

As well as R and python (2.7+) which most of the code is written in, you will need the packages spindrift (https://github.com/mathii/spindrift) and pyEigenstrat (https://github.com/mathii/pyEigenstrat).

- Iain Mathieson 15Oct15
