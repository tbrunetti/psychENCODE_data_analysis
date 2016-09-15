import os
import re
import subprocess
from chunkypipes.components import *

class Pipeline(BasePipeline):
	def dependencies(self):
		return []

	def description(self):
		return{'Run IDR for individual replication consistency'}

	def configure(self):
		return {
			'samtools':{
				'path':'Full path to samtools exectuable'
			},
			'bedtools':{
				'path':'Full path to bedtools bamToBed exectuable'
			},
			'idr':{
				'path':'Full path to IDR scripts, full path must end in /'
			},
			'macs2':{
				'path':'Full path to macs2'
			}
		}

	def add_pipeline_args(self, parser):
		parser.add_argument('-inputBAM', required=True, help="Full path to BAM file to be analyzed")
		parser.add_argument('-narrowPeaks', required=True, help='Full path to narrowPeaks file from MACS2 of inputBAM')
		parser.add_argument('-outDir', default=os.getcwd(), help='Full path to output directory, ending in /')


	def run_pipeline(self, pipeline_args, pipeline_config):
		samtools_sort=Software('samtools_sort', pipeline_config['samtools']['path']+' sort')
		samtools_view=Software('samtools_view', pipeline_config['samtools']['path']+' view')
		bedtools_bamToBed=Software('bedtools_bamToBed', pipeline_config['bedtools']['path'])
		macs2_callpeak=Software('macs2_callpeak', pipeline_config['macs2']['path']+' callpeak')

		# regex to extract all characters after the last / in input file
		bamNameRegex=re.compile(r'([^/]+$)')
		bamName=re.search(bamNameRegex, str(pipeline_args['inputBAM'])).group(1)
		
		bedConversion=open(pipeline_args['outDir']+str(bamName[:-4])+'_bed_converted.bedpe', 'w')
		tagAlignFile=open(pipeline_args['outDir']+str(bamName[:-4])+'.bedpe.tagAlign', 'w')
		
		# filters bams prior to sorting for BEDPE conversion for IDR use
		samtools_view.run(
			Parameter('-uf', '0x2'),
			Parameter('-F', '1548'),
			Parameter('-q', '30'),
			Parameter(pipeline_args['inputBAM']),
			Pipe(
				samtools_sort.pipe(
					Parameter('-n'),
					Parameter('-'),
					Parameter(pipeline_args['outDir']+str(bamName[:-4])+'_sortedByName')
				)
			)
		)
		
		# converts BAM to BEDPE for creation of pseudoreplicates
		bedtools_bamToBed.run(
					Parameter('-i', pipeline_args['outDir']+str(bamName[:-4])+'_sortedByName.bam'),
					Parameter('-bedpe'),
					Redirect(stream=Redirect.STDOUT, dest=pipeline_args['outDir']+str(bamName[:-4])+'_bed_converted.bedpe')
		)



		# makes bedpe into a "bed" where left most pair-end tag is start, and right pair-end tag is end
		with open(pipeline_args['outDir']+str(bamName[:-4])+'_bed_converted.bedpe') as convertToBed:
			for line in convertToBed:
				if line=='\n':
					continue;
				else:
					chrpos1, start1, end1, chrpos2, start2, end2, name, score, strand1, strand2=line.split('\t')
					bedformat=[chrpos1, start1, end2, name, score, strand1, strand2]
					tagAlignFile.write('\t'.join(bedformat))
		
		subprocess.call(['gzip', '-c', pipeline_args['outDir']+str(bamName[:-4])+'.bedpe.tagAlign'], stdout=open(pipeline_args['outDir']+str(bamName[:-4])+'.bedpe.tagAlign.gz', 'w'))
		# create two pseudoreplicates
		open_tagAlign = subprocess.Popen(['zcat', pipeline_args['outDir']+str(bamName[:-4])+'.bedpe.tagAlign.gz'], stdout=subprocess.PIPE)
		get_num_reads = subprocess.Popen(['wc', '-l'], stdin = open_tagAlign.stdout, stdout=subprocess.PIPE)
		

		# total reads in a file
		total_reads = get_num_reads.stdout.read()
		# number of reads in each pseudo replicate
		reads_per_pseudo_rep = int(total_reads)/2

		# shuffle file for randomness and split into two pseudoreplicates
		open_tagAlign = subprocess.Popen(['zcat', pipeline_args['outDir']+str(bamName[:-4])+'.bedpe.tagAlign.gz'], stdout=subprocess.PIPE)
		shuffle_reads = subprocess.Popen(['shuf'], stdin=open_tagAlign.stdout, stdout=subprocess.PIPE)
		split_reads = subprocess.Popen(['split', '-d', '-l', str(reads_per_pseudo_rep), '-', str(pipeline_args['outDir']+str(bamName[:-4]))], stdin=shuffle_reads.stdout)
		# waits for subprocess to finish running
		split_reads.communicate()
		
		# gzip and rename split reads pseudo replicate tagAlign files
		subprocess.call(['gzip', str(pipeline_args['outDir']+str(bamName[:-4]))+'00'])
		subprocess.call(['gzip', str(pipeline_args['outDir']+str(bamName[:-4]))+'01'])
		subprocess.call(['mv', str(pipeline_args['outDir']+str(bamName[:-4]))+'00'+'.gz', str(pipeline_args['outDir']+str(bamName[:-4]))+'.pr1.bedpe.tagAlign.gz'])
		subprocess.call(['mv', str(pipeline_args['outDir']+str(bamName[:-4]))+'01'+'.gz', str(pipeline_args['outDir']+str(bamName[:-4]))+'.pr2.bedpe.tagAlign.gz'])

		
		# call peaks on pseudoreplicates under same conditions as original file
		macs2_callpeak.run(
				Parameter('-t', str(pipeline_args['outDir']+str(bamName[:-4]))+'.pr1.bedpe.tagAlign.gz'),
				Parameter('-f', 'BEDPE'),
				Parameter('-g', 'hs'),
				Parameter('-n', str(pipeline_args['outDir']+str(bamName[:-4]))+'.pr1.bedpe.tagAlign' + '_regular_peak_calls'),
				Parameter('--nomodel'),
				Parameter('--extsize', '200'),
				Parameter('--shift', '-100'),
				Parameter('-B', '--SPMR'), # Generates pileup tracks, bedgraph, fragment pileup per million reads
				Parameter('--call-summits'),
				Parameter('--keep-dup', 'all')
			)


		macs2_callpeak.run(
				Parameter('-t', str(pipeline_args['outDir']+str(bamName[:-4]))+'.pr2.bedpe.tagAlign.gz'),
				Parameter('-f', 'BEDPE'),
				Parameter('-g', 'hs'),
				Parameter('-n', str(pipeline_args['outDir']+str(bamName[:-4]))+'.pr2.bedpe.tagAlign' + '_regular_peak_calls'),
				Parameter('--nomodel'),
				Parameter('--extsize', '200'),
				Parameter('--shift', '-100'),
				Parameter('-B', '--SPMR'), # Generates pileup tracks, bedgraph, fragment pileup per million reads
				Parameter('--call-summits'),
				Parameter('--keep-dup', 'all')
			)

		os.chdir(pipeline_config['idr']['path'])
		# pairwise comparison of peaks calls: original vs pseudorep1, orginal vs pseudorep2, pseudorep1 vs pseudorep2
		subprocess.call(['Rscript', pipeline_config['idr']['path']+'batch-consistency-analysis.r', pipeline_args['narrowPeaks'], str(pipeline_args['outDir']+str(bamName[:-4]))+'.pr1.bedpe.tagAlign' + '_regular_peak_calls_peaks.narrowPeak', '-1', str(pipeline_args['outDir']+str(bamName[:-4]))+'_bedpe_original.vs.pr1', '0', 'F', 'p.value'])
		subprocess.call(['Rscript', pipeline_config['idr']['path']+'batch-consistency-analysis.r', pipeline_args['narrowPeaks'], str(pipeline_args['outDir']+str(bamName[:-4]))+'.pr2.bedpe.tagAlign' + '_regular_peak_calls_peaks.narrowPeak', '-1', str(pipeline_args['outDir']+str(bamName[:-4]))+'_bedpe_original.vs.pr2', '0', 'F', 'p.value'])
		subprocess.call(['Rscript', pipeline_config['idr']['path']+'batch-consistency-analysis.r', str(pipeline_args['outDir']+str(bamName[:-4]))+'.pr1.bedpe.tagAlign' + '_regular_peak_calls_peaks.narrowPeak', str(pipeline_args['outDir']+str(bamName[:-4]))+'.pr2.bedpe.tagAlign' + '_regular_peak_calls_peaks.narrowPeak', '-1', str(pipeline_args['outDir']+str(bamName[:-4]))+'_bedpe_pr1.vs.pr2', '0', 'F', 'p.value'])

		# plot all pairwise comparison using IDR plot
		subprocess.call(['Rscript', pipeline_config['idr']['path']+'batch-consistency-plot.r', '3', pipeline_args['outDir']+str(bamName[:-4]), str(pipeline_args['outDir']+str(bamName[:-4]))+'_bedpe_original.vs.pr1', str(pipeline_args['outDir']+str(bamName[:-4]))+'_bedpe_original.vs.pr2', str(pipeline_args['outDir']+str(bamName[:-4]))+'_bedpe_pr1.vs.pr2'])
		
		# extract peaks with IDR < 0.05? 0.01?
		passed_IDR_origpr1=[]
		passed_IDR_origpr2=[]
		passed_IDR_pr1Vspr2=[]
		with open(pipeline_args['outDir']+str(bamName[:-4])+'_bedpe_original.vs.pr1-overlapped-peaks.txt') as originalpr1:
			header=next(originalpr1)
			for line in originalpr1:
				line=line.split(' ')
				if float(line[10].rstrip('\n')) <= 0.01:
					passed_IDR_origpr1.append(line)
				else:
					continue;

		with open(pipeline_args['outDir']+str(bamName[:-4])+'_bedpe_original.vs.pr2-overlapped-peaks.txt') as originalpr2:
			header=next(originalpr2)
			for line in originalpr2:
				line=line.split(' ')
				if float(line[10].rstrip('\n')) <= 0.01:
					passed_IDR_origpr2.append(line)
				else:
					continue;

		with open(pipeline_args['outDir']+str(bamName[:-4])+'_bedpe_pr1.vs.pr2-overlapped-peaks.txt') as pr1Vspr2:
			header=next(pr1Vspr2)
			for line in pr1Vspr2:
				line=line.split(' ')
				if float(line[10].rstrip('\n')) <= 0.01:
					passed_IDR_pr1Vspr2.append(line)
				else:
					continue;

		f1=open(pipeline_args['outDir']+str(bamName[:-4])+'_bedpe_original.vs.pr1_peaks_passed_IDR.txt', 'w')
		f2=open(pipeline_args['outDir']+str(bamName[:-4])+'_bedpe_original.vs.pr2_peaks_passed_IDR.txt', 'w')
		f3=open(pipeline_args['outDir']+str(bamName[:-4])+'_bedpe_pr1.vs.pr2_peaks_passed_IDR.txt', 'w')
		f4=open(pipeline_args['outDir']+str(bamName[:-4])+'_bedpe_stats_peaks_passed_IDR.txt', 'w')
		
		for x in range(0, len(passed_IDR_origpr1)):
			f1.write('\t'.join(passed_IDR_origpr1[x]))

		for n in range(0, len(passed_IDR_origpr2)):
			f2.write('\t'.join(passed_IDR_origpr2[n]))

		for z in range(0, len(passed_IDR_pr1Vspr2)):
			f3.write('\t'.join(passed_IDR_pr1Vspr2[z]))

				
		numOrigPR1=subprocess.check_output(['wc', '-l', pipeline_args['outDir']+str(bamName[:-4])+'_bedpe_original.vs.pr1_peaks_passed_IDR.txt'])
		numOrigPR2=subprocess.check_output(['wc', '-l', pipeline_args['outDir']+str(bamName[:-4])+'_bedpe_original.vs.pr2_peaks_passed_IDR.txt'])
		numPR1vsPR2=subprocess.check_output(['wc', '-l', pipeline_args['outDir']+str(bamName[:-4])+'_bedpe_pr1.vs.pr2_peaks_passed_IDR.txt'])


		totalPeaks_orig = subprocess.check_output(['wc', '-l', pipeline_args['narrowPeaks']])
		pseudoreplicates1 = subprocess.check_output(['wc', '-l', str(pipeline_args['outDir']+str(bamName[:-4]))+'.pr1.bedpe.tagAlign' + '_regular_peak_calls_peaks.narrowPeak'])
		pseudoreplicates2 = subprocess.check_output(['wc', '-l', str(pipeline_args['outDir']+str(bamName[:-4]))+'.pr2.bedpe.tagAlign' + '_regular_peak_calls_peaks.narrowPeak'])
		
		# writes results to file
		f4.write('sample_name'+'\t'+str(bamName)+'\n'+'\n')
		
		f4.write('-----------------Total_Peaks_Across_All_Files-----------------'+'\n')
		f4.write('total_original_peaks'+'\t'+str(totalPeaks_orig)+'\n')
		f4.write('total_pseudoRep1_peaks'+'\t'+str(pseudoreplicates1)+'\n')
		f4.write('total_pseudoRep2_peaks'+'\t'+str(pseudoreplicates2)+'\n'+'\n')

		f4.write('-----------------Total_Peaks_Passed_IDR_Threshold_Across_All_Files-----------------'+'\n')
		f4.write('original_vs_pseudoRep1'+'\t'+str(numOrigPR1)+'\n')
		f4.write('original_vs_pseudoRep2'+'\t'+str(numOrigPR2)+'\n')
		f4.write('pseudoRep1_vs_pseudoRep2'+'\t'+str(numPR1vsPR2))