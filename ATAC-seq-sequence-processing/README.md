#ATAC-seq-analysis-pipeline-version2.py and ATAC-seq-analysis-pipeline-version3-bedpe.py
------------------------------------------------------------------------------------------
A pipeline to process ATAC-seq data and generate peak files of chromatin accessibility.  In version2, peaks called are processed using single-end reads despite inputing paired-end data.  Version3-bedpe, generates peaks using paired-end data information.

###Overview
-----------
![Alt text](https://github.com/tbrunetti/psychENCODE_data_analysis/blob/master/ATAC-seq-sequence-processing/ATAC-seq-pipeline-overview.jpg)

###Software Requirements
-------------------------
* Python minimum version requirement 2.7.6 
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

###User Generated/Provided File Requirements
----------------------------------------------
* fasta of reference genome indexed by bwa
* bed file of blacklisted genomic regions
  * This file can be downloaded from mod/mouse/humanENCODE
  * https://sites.google.com/site/anshulkundaje/projects/blacklists
* genome size file (i.e. chrID followed by size_of_chromosome)
  * This file can be downloaded from UCSC

###Installation and Configuration
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
Post-installation, the file ATAC-seq-analysis-pipeline-version3-bedpe.pyc is created.  This is the configuration file for the installed pipeline.  In order to configure this file, run the following command:
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
The user will notice that once the information is recorded, chunkypipes will state "Configuration file successfully written" to notify the user the configuration has been successfully modified and  saved. This file only needs to be configured once, unless the paths to these files has changed or if a different number of threads is to be used.

###Running the Pipeline
-------------------------
Both pipelines are equipped with a -h or --help flag for details on default parameters and user options.  An example of how to use this is shown below:
```
chunky run ATAC-seq-analysis-pipeline-version3-bedpe.py -h
```
The minimum required arguments is the --reads flag and the --output flag. Reads needs to be in fastq format with the name of the forward pairs followed by a colon followed by the name of the respective reverse pairs.  The output should be the full path to the desired user generated output directory.
```
chunky run ATAC-seq-analysis-pipeline-version3-bedpe.py --reads forwardReads.fastq:reverseReads.fastq --output ~/my_results
```
__***Optional Arguments***__
* **--lib**   *Specifiy the name of sample; default is wall-time as sample name*
* **--forward-adapter**   *specify adapter sequence on forward strand; default is none*
* **--reverse-adapter**   *specify adapter sequence on reverse strand; default is none*

__***Default Software Parameters***__
* cutadapt
  * minimum read length of 5 bp
  * minimum quality score of 30 
* bwa
  * maximum alignments to output is 1
  * maximum insert size for proper mapping is 2000
* samtools
  * *For uniquely mapped reads:*
    * -F 256 (not primary alignment)
    * skip alignments with MAPQ smaller than 10
  * *For removing unmapped reads:*
    * -F 12 (read unmapped/mate unmapped)
  * *For converting BAM to BEDPE*
    * -uf 0x2 (properly paired reads only)
    * -F 1548 (read/mate unmapped, fails quality checks, read is PCR/optical duplicate)
    * skip alignments with score less than 30
* bedtools
  * *intersect (for use in blacklisted genomic regions)*
    * -v
    * -f 0.5 (50% minimum overlap required)
* Picard
  * REMOVE_DUPLICATES=true
  * VALIDATION_STRINGENCY=LENIENT
* MACS2   
MACS2 is called a total of 4 times.  Twice for BED and twice for BEDPE.  Note this is NOT the same BEDPE adopted by most genome centers.  This BEDPE is specially formatted for MACS2 specifications.

  * *narrow peaks (default q-value=0.01)*
    * --nomodel
    * --extsize 200
    * --shift -100
    * -B --SPMR (bedgraph pileup tracks)
    * --call-summits
    * --keep-dup all
  * *broad peaks*
    * -q 0.05
    * --nomodel
    * --extsize 200
    * --shift -100
    * --broad
    * --keep-dup all

###Pipeline Output
------------------
__***sample_name.<int>.bam***__ is the most unprocessed bam file with no altercations or cleaning  
__***sample.name.processed.bam***__ is the most processed bam file and is the bam ultimately used for all end analysis  
__***sample.name.unique.bam***__ is the bam where everything is processed down to the unique reads except unmapped reads, chrM, and blacklisted regions  have not yet been removed.  

The following bam files are the binned reads and each has their own set of flagstats:
* sample_name.chrM
* sample_name.nucleosome_free
* sample_name.mononucleosome
* sample_name.dinucleosome
* sample_name.trinucleosome

__***sample_name.unshifted***__ are reads that come from the final processed bam file but have not yet been shifted to correct for transposase cutting  

__***fastqc directory***__ directory containing all the FastQC results  
__***logs directory***__ directory that contains quality control information and picard statistics

