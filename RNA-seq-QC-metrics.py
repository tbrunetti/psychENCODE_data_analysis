import os
from chunkypipes.components import *
import subprocess
import datetime

#not to be used for raw sequences!  Pipeline QCs aligned bam files only!
#note picard and samtools should be installed prior to running pipeline
class Pipeline(BasePipeline):

	def dependencies(self):
		#use pip compatible install name
		#RSeQC requires dependencies: bx-python and pysam but 
		#pip install RSeQC should automatically install these if not already installed
		return ['RSeQC']

	def description(self):
		return 'QC files for all RNA-seq bams as input'

	def configure(self):
		return{
			'general':{
				'reference_genome':"Full path to reference genome in fasta format",
			},

			'samtools':{
				'path':"Full path to samtools executable",
				'output':'Full path to create output directory'
			},

			'picard':{
				'path':"Full path to picard jar file, i.e. java -jar $/Tools/picardTools/picard.jar",
				'ref_flat':"Full path to gene annotations in refFlat form, REQUIRED",
				'rRNA_interval_list':"Full path to rRNA interval list",
				'output':"Full path to create output directory of Picard results"
			},

			'RSeQC':{
				'pathTotin':"Full path to RSeQC tin.py script",
				'pathToBamstats':'Full pathe to RSeQC bam_stat.py script',
				'refSeqGenes':"Full path to reference gene model in BED format (12-column BED file)",
				'pathToAlignStats':"Full path to create output directory of alignment statistics, should be separate from samtools output"
			}
		}

	def add_pipeline_args(self, parser):
		#arguments for samtools
		parser.add_argument('-q', default='20', help="[INT] min MAPQ score to be considered unique read, default=20")
		parser.add_argument('-b', required=True, help="Full path to directory of aligned bams")
		
		#arguments for picard CollectRnaSeqMetrics
		parser.add_argument('-STRAND', required=True, help='options: NONE, FIRST_READ_TRANSCRIPTION_STRAND, or SECOND_READ_TRANSCRIPTION_STRAND')
		parser.add_argument('-sort', default='FALSE', help='[boolean] True or False, default=False')

	def run_pipeline(self, pipeline_args, pipeline_config):
		
		#this function checks if a directory name already exists, if yes, then it
		#continues to check for exisitng names recursively until a name has not yet
		#been used
		def making_directories(directoryNumber):
			#name of directories to make
			pathToSamOut=pipeline_config['samtools']['output']+'RNAseq.bam.output'+'_dir_'+str(directoryNumber)+'/'
			pathToPicardOut=pipeline_config['picard']['output']+'picard.RNAstats'+'_dir_'+str(directoryNumber)+'/'
			pathToRseqcOut=pipeline_config['RSeQC']['pathToAlignStats']+'RSeQC.bamAlignStats'+'_dir_'+str(directoryNumber)+'/'
			#check if directory already exist before creating it
			if os.path.exists(pathToSamOut)==False:
				os.mkdir(pathToSamOut)
				os.mkdir(pathToPicardOut)
				os.mkdir(pathToRseqcOut)
				return directoryNumber, pathToSamOut, pathToPicardOut, pathToRseqcOut
			else:
				directoryNumber+=1
				return making_directories(directoryNumber)
		
		#always start at directory number 0	and call making_directories for output
		directoryNumber, pathToSamOut, pathToPicardOut, pathToRseqcOut=making_directories(0)
		
		samtools=Software('samtools', pipeline_config['samtools']['path'])
		picard=Software('picard', pipeline_config['picard']['path'])
		rseqcBamStats=Software('RSeQC', pipeline_config['RSeQC']['pathToBamstats'])
		rseqcTin=Software('RSeQC', pipeline_config['RSeQC']['pathTotin'])
		
		#-------filter all bam files for uniquely mapped reads only------------
		
		for bam in os.listdir(pipeline_args['b']):
			#must use this check because .bai files will exist in same directory
			if bam[-4:]=='.bam':
				print "extracting unique reads from "+str(bam)
				try:
					samtools.run(
						#runs view method and adds header to new bam file
						Parameter('view', '-h'),
						Parameter('-q', pipeline_args['q']),
						Parameter('-b', pipeline_args['b']+str(bam)),
						#samflag 256=not primary alignment
						Parameter('-F', '256'),
						Redirect(stream=Redirect.STDOUT, dest=pathToSamOut+'unique.'+str(bam))
					)
				
					#create new bam index corresponding to the unique reads bam file
					samtools.run(
						Parameter('index'),
						#output of unique bams becomes input into index call
						Parameter(pathToSamOut+'unique.'+str(bam)),
						Redirect(stream=Redirect.STDOUT, dest=pathToSamOut+'unique.'+str(bam)+'.bai')
					)

				except OSError:
					print 'OSError: refer to Python OSError exception for causes of failure'
					sys.exit()

				except:
					print 'Error: ' + str(bam) + ' could not be processed for extraction and .bai was not created'	
		
		#----------------run picard collectRnaSeqMetrics-----------------------
		
		#creates a cumulative text file of RNA metrics for all bams analyzed
		filename='RNA-seq-sample-metrics_'+str(directoryNumber)+'.txt'
		sampleMetrics=open(filename, 'w')
		for bam in os.listdir(pipeline_args['b']):
			if bam[-4:]=='.bam':
				print "Running Picard CollectRnaSeqMetrics on "+str(bam)
				picard.run(
					Parameter('CollectRnaSeqMetrics'),
					Parameter('REF_FLAT='+pipeline_config['picard']['ref_flat']),
					Parameter('RIBOSOMAL_INTERVALS='+pipeline_config['picard']['rRNA_interval_list']),
					Parameter('STRAND='+pipeline_args['STRAND']),
					Parameter('INPUT='+pathToSamOut+'unique.'+str(bam)),
					Parameter('OUTPUT='+pathToPicardOut+'unique.'+str(bam[:-4])+'.rnaSeqMetrics.txt'),
					Parameter('ASSUME_SORTED='+pipeline_args['sort'])
					)
				print "Finished running picard, concatenating result to unified file"
				#unify all sample results into one file
				f=open(pathToPicardOut+'unique.'+str(bam[:-4])+'.rnaSeqMetrics.txt')
				line=f.readlines()
				header=line[6].split('\t')
				metricsInfo=line[7].split('\t')
				metricsInfo[22]=str(bam[:-4])
				with open(filename, 'a+') as file:
					if os.stat(filename).st_size==0:
						for labels in range(0, len(header)-1):
							file.write(str(header[labels])+'\t')
						file.write(header[len(header)-1])
					
					for col in range(0, len(metricsInfo)-1):
						file.write(str(metricsInfo[col])+'\t')
					file.write(str(metricsInfo[len(metricsInfo)-1]))


		
		#---------run RSeQC to get general BAM alignment statistics--------

		#calculates and outputs alignment statistics from BAM files
		for bam in os.listdir(pathToSamOut):
			if bam[-4:]=='.bam':
				print "Calculating alignment statistic for "+str(bam)
				rseqcBamStats.run(
					Parameter('-i', pathToSamOut+str(bam)),
					Parameter('-q', pipeline_args['q']),
					Redirect(stream=Redirect.STDOUT, dest=pathToRseqcOut+str(bam)+'.alignStats')
				)


		#----------------run RSeQC transcript integrity number (TIN)----------

		#calculates and outputs TIN scores fore each transcript type in samples and an
		#estimated RIN score for the entire sample
		#COMPUTATIONAL BOTTLE NECK!!!!!! Comment run line below if TIN and RIN scores not needed
		print "Running bams for transcript integrity number calculation"
		rseqcTin.run(
			Parameter('-i', pathToSamOut),
			Parameter('-r', pipeline_config['RSeQC']['refSeqGenes'])
			)