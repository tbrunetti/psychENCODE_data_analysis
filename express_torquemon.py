import os

# specific function for eXpress RNA-seq quantifier
# creates new output folders for each sample
def format_outputDir_express():
	os.chdir(pathToQuantDir)
	for directories in os.listdir(pathToRNAseqDirectories):
		os.mkdir('express-'+str(directories))

# generation of tqm file for Beagle job submission
# configured for eXpress with input from STAR trancriptome
def generate_torquemon_file():
	torquemonFile = open('RNAseq-quant.tqm', 'w')
	os.chdir(pathToRNAseqDirectories)
	for directories in os.listdir(pathToRNAseqDirectories):
		os.chdir(directories)
		cwd = os.getcwd()
		inputBam = str(cwd)+'/'+str(directories)+'_0.Aligned.toTranscriptome.out.bam'
		outputDir = pathToQuantDir+'express-'+str(directories)
		torquemonFile.write('qsub -v inputBam='+str(inputBam)+',outputDir='+str(outputDir)+' '+str(pathToPBS)+'\t'+str(outputDir)+'\n')
		os.chdir(pathToRNAseqDirectories)

if __name__=='__main__':
	pathToQuantDir = '/lustre/beagle/tbrunetti/RNA_quantification/quantification_express/'
	pathToRNAseqDirectories = '/lustre/beagle2/djf604/projects/BrainGVEX_local/storage/'
	pathToPBS = '/lustre/beagle/tbrunetti/RNA_quantification/express-quant.pbs'
	
	#format_outputDir_express();
	generate_torquemon_file();