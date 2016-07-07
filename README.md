# psychENCODE_data_analysis

The scripts and programs listed in this repository have been designed for use with the psychENCODE projects.  Most of these scripts can be adapted for use on general projects pertaining to RNA-seq QC and analysis.

##RNA-seq-QC-metrics.py
-----------------------
This is a pipline that is designed to use various metrics to obtain quality control data from RNA-seq related experiments.  It does not yet having visualization functionality for all experiments run, but there is a plan to implement the output of the data generated from this pipeline into an R script for statistics and visualization.

####Software Requirements
-------------------------
* chunkypipes (http://chunky-pipes.readthedocs.io/en/stable/getting_started.html)
* samtools (http://samtools.sourceforge.net)
* Picard (https://broadinstitute.github.io/picard)
* RSeQC (http://rseqc.sourceforge.net)
 * Upon pipeline configuration, this requirement can be satisfied automatically


####User generated File Requirements
------------------------------------
In addtion to the software requirements stated above, the following files must be generated and provided by the user:
* reference genome in a single FASTA file, as well as corresponding dictionary and .fai files in the same directory
* rRNA coordiates in interval list format
* genome annotations in refFlat format
* reference gene model in BED-12 format

####Installation and Configuration
----------------------------------
Assuming chunkypipes has been installed correctly, download RNA-seq-QC-metrics.py and run the following:

```
chunky install RNA-seq-QC-metrics.py
```
Upon successful installation, it should then ask the user whether the RSeQC dependecy should be intalled. If this dependency has already been installed the user can skip this step by selecting 'n' for no.

```
Pipeline RNA-seq-QC-metrics.py successfully installed.

Attempting to install the following dependencies:
RSeQC

Proceed with dependency installation? [y/n] 
```

Now the pipleline is ready for configuration.  Run the following command to configure the pipeline:
```
chunky configure RNA-seq-QC-metrics.py
```
This will then prompt the user to input the paths to a number of tools and files.
