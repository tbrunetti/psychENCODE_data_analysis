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
* bedtools; must be version >=2.25.0 (http://bedtools.readthedocs.io/en/latest/)
* macs2 (https://pypi.python.org/pypi/MACS2)

