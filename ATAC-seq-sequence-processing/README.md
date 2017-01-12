##ATAC-seq-analysis-pipeline-version2.py and ATAC-seq-analysis-pipeline-version3-bedpe.py
------------------------------------------------------------------------------------------
A pipeline to process ATAC-seq data and generate peak files of chromatin accessibility.  In version2, peaks called are processed using single-end reads despite inputing paired-end data.  Version3-bedpe, generates peaks using paired-end data information.

####Software Requirements
-------------------------
* chunkypipes (http://chunky-pipes.readthedocs.io/en/stable/getting_started.html)
* cutadapt (http://cutadapt.readthedocs.io/en/stable/installation.html)
* fastqc (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* bwa (http://bio-bwa.sourceforge.net/)
* novosort (http://www.novocraft.com/products/novosort/)
* picard (https://broadinstitute.github.io/picard/)
* samtools (http://samtools.sourceforge.net)
* bedtools (http://bedtools.readthedocs.io/en/latest/)
  * must be version >=2.25.0
* macs2 (https://pypi.python.org/pypi/MACS2)

####User Generated/Provided File Requirements
----------------------------------------------
* fasta of reference genome indexed by bwa
* bed file of blacklisted genomic regions
  * This file can be downloaded from mod/mouse/humanENCODE
  * https://sites.google.com/site/anshulkundaje/projects/blacklists
* genome size file (i.e. chrID followed by size_of_chromosome)
  * This file can be downloaded from UCSC

####Installation and Configuration
----------------------------------
Either download the pipleline version that would like to be used or clone the repository into your own directory. 
```
mkdir ~/my_project
cd ~/my_project
git clone https://github.com/tbrunetti/psychENCODE_data_analysis.git
cd ATAC-seq-sequence-processing
```
Assuming chunkypipes has been installed correctly, run the following:

```
chunky install ATAC-seq-analysis-pipeline-version3-bedpe.py

```
If the installation was successful the following message will appear:
```
Pipeline ATAC-seq-analysis-pipeline-version3-bedpe.py successfully installed.
```
Additionally, you should notice post installation the creation of the file ATAC-seq-analysis-pipeline-version3-bedpe.pyc.  This is the configuration file for the installed pipeline.  In order to configure the pipeline, run the following command:
```
chunky configure ATAC-seq-analysis-pipeline-version3-bedpe.py
```
This will prompt the user for the following infomration:  NOTE! It is critical that the FULL FILE PATH is written out at each prompt.
```
Full path to novosort []: 
Number of threads to use for Novosort []: 
Full path to bwa executable []: 
Number of threads to use for bwa aln []: 
Directory of the bwa reference index [Ex. /path/to/bwa/index/genome.fa] []: 
Full path to samtools []: 
Full path to FastQC []: 
Full path to MACS2 []: 
Full path to cutadapt executable []: 
Full path to Picard [Ex. java -jar /path/to/picard.jar] []: 
Full path to bedtools >=2.25.0 []: 
Full path to a genome sizes file []: 
Full path to the BED of blacklisted genomic regions []: 
Configuration file successfully written.
```
The user will notice that once the information is recorded, chunkypipes will state "Configuration file successfully written" to notify the user the configuration has been successfully modified and  saved.
