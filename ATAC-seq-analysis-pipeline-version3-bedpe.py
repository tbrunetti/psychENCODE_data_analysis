import os
import subprocess
import re
import time
import pysam
from chunkypipes.components import Software, Parameter, Redirect, Pipe, BasePipeline

#credited to Dominic Fitzgerald (djf604)
#modified by Tonya Brunetti (tbrunetti): peak calling method, binning of reads, and added some metric collection points
READ1 = 0
FIRST_CHAR = 0

MINUS_STRAND_SHIFT = -5
PLUS_STRAND_SHIFT = 4

STERIC_HINDRANCE_CUTOFF = 38

class Pipeline(BasePipeline):
	def description(self):
		return ['Pipeline for analyzing ATAC-seq data']

	def configure(self):
		return{
			'cutadapt':{
				'path': 'Full path to cutadapt executable'
			},
			'bwa':{
				'path': 'Full path to bwa executable',
				'threads': 'Number of threads to use for bwa aln',
				'index-dir': 'Directory of the bwa reference index [Ex. /path/to/bwa/index/genome.fa]'
			},
			'fastqc':{
				'path': 'Full path to FastQC'
			},
			'samtools':{
				'path': 'Full path to samtools'
			},
			'novosort':{
				'path': 'Full path to novosort',
				'threads': 'Number of threads to use for Novosort'
			},
			'picard':{
				'path': 'Full path to Picard [Ex. java -jar /path/to/picard.jar]'
			},
			'bedtools':{
				'path': 'Full path to bedtools >=2.25.0',
				'blacklist-bed': 'Full path to the BED of blacklisted genomic regions',
				'genome-sizes': 'Full path to a genome sizes file'
			},
			'macs2':{
				'path': 'Full path to MACS2'
			}
		}

	def add_pipeline_args(self, parser):
		parser.add_argument('--reads', required=True, help='read1:read2', action='append')
		parser.add_argument('--output', required=True)
		parser.add_argument('--lib', default=str(time.time()))
		parser.add_argument('--step', default=0)
		parser.add_argument('--forward-adapter', default='ZZZ')
		parser.add_argument('--reverse-adapter', default='ZZZ')
		return parser

	@staticmethod
	def count_gzipped_lines(filepath):
		zcat = subprocess.Popen(['zcat', filepath], stdout=subprocess.PIPE)
		num_lines = subprocess.check_output(['wc', '-l'], stdin=zcat.stdout)
		return num_lines.strip()

	@staticmethod
	def shift_reads(input_bed_filepath, output_bed_filepath, genome_sizes_filepath,
					log_filepath, minus_strand_shift, plus_strand_shift):
		def to_log(msg, log_filepath):
			with open(log_filepath, 'a') as log_file:
				log_file.write(msg + '\n')

		def within_genome_boundary(chrom, check_pos, genome_sizes_dict):
			if check_pos < 1:
				return False
			if check_pos > genome_sizes_dict[chrom]:
				return False
			return True

		genome_sizes={}
		with open(genome_sizes_filepath) as genome_sizes_file:
			for line in genome_sizes_file:
				chrom, size = line.strip().split('\t')
				genome_sizes[chrom] = int(size)

		output_bed = open(output_bed_filepath, 'w')
		with open(input_bed_filepath) as input_bed:
			for line in input_bed:
				record=line.strip('\n').split('\t')
				chrom, start, end, strand = record[0], int(record[1]), int(record[2]), record[5]

				if strand == '+':
					new_start = start + plus_strand_shift
					new_end = end + plus_strand_shift
				elif strand == '-':
					new_start = start + minus_strand_shift
					new_end = end + minus_strand_shift
				else:
					to_log('Mal-formatted line in file, skipping:', log_filepath)
					to_log(line.strip(), log_filepath)
					continue

				if not within_genome_boundary(chrom, new_start, genome_sizes):
					to_log('Shift end site beyond chromosome boundary:', log_filepath)
					to_log(line.strip(), log_filepath)
					continue

				output_bed.write('\t'.join([chrom, str(new_start), str(new_end)] + record[3:]) + '\n')

			output_bed.flush()
			output_bed.close()

	@staticmethod
	# slightly different shift method for bedpe format
	def shift_reads_bedpe(input_bed_filepath, output_bed_filepath, genome_sizes_filepath,
					log_filepath, minus_strand_shift, plus_strand_shift):

		# next two def are the same as above shift_reads method
		def to_log(msg, log_filepath):
			with open(log_filepath, 'a') as log_file:
				log_file.write(msg + '\n')

		def within_genome_boundary(chrom, check_pos, genome_sizes_dict):
			if check_pos < 1:
				return False
			if check_pos > genome_sizes_dict[chrom]:
				return False
			return True

		genome_sizes={}
		with open(genome_sizes_filepath) as genome_sizes_file:
			for line in genome_sizes_file:
				chrom, size = line.strip().split('\t')
				genome_sizes[chrom] = int(size)

		output_bed = open(output_bed_filepath, 'w')
		with open(input_bed_filepath) as input_bed:
			for line in input_bed:
				record=line.strip('\n').split('\t')
				if len(record)<7:
					to_log('Mal-formatted line in file, skipping:', log_filepath)
					to_log(line.strip(), log_filepath)
					continue;

				else:
					chrom, start, end, name, score, strand1, strand2 = record[0], int(record[1]), int(record[2]), record[3], int(record[4]), record[5], record[6].rstrip('\n')

					if strand1 == '+' and strand2 == '-':
						new_start = start + plus_strand_shift
						new_end = end + minus_strand_shift

				
					elif strand1 == '-' and strand2 == '+':
						new_start = start + minus_strand_shift
						new_end = end + plus_strand_shift
				
					else:
						to_log('Mal-formatted line in file, skipping:', log_filepath)
						to_log(line.strip(), log_filepath)
						continue

					if not within_genome_boundary(chrom, new_start, genome_sizes):
						to_log('Shift end site beyond chromosome boundary:', log_filepath)
						to_log(line.strip(), log_filepath)
						continue

					output_bed.write('\t'.join([chrom, str(new_start), str(new_end)] + record[3:]) + '\n')

			output_bed.flush()
			output_bed.close()





	def run_pipeline(self, pipeline_args, pipeline_config):
		# Instantiate variable from argparse
		read_pairs = pipeline_args['reads']
		output_dir = os.path.abspath(pipeline_args['output'])
		logs_dir = os.path.join(output_dir, 'logs')
		lib_prefix = pipeline_args['lib']
		step = int(pipeline_args['step'])
		forward_adapter = pipeline_args['forward_adapter']
		reverse_adapter = pipeline_args['reverse_adapter']

		# Create output, tmp, and logs directories
		tmp_dir = os.path.join(output_dir, 'tmp')
		subprocess.call(['mkdir', '-p', output_dir, tmp_dir, logs_dir])

		#Keep list of items to delete
		staging_delete = [tmp_dir]
		bwa_bam_outs = []
		qc_data = {
			'total_raw_reads_counts': [],
			'trimmed_reads_counts': [],
			'num_reads_mapped': [],
			'num_read_removed_steric_hinderence': '0',
			'percent_duplicate_reads': '0',
			'num_unique_reads_mapped': [], #implemented
			'num_mtDNA_reads_mapped': [],
			'percent_mtDNA_reads_mapped': '0' ,
			'num_reads_mapped_after_filtering': '-1', #TODO This isn't implemented
			'num_peaks_called': '-1',
			#TODO Get number of peaks in annotation sites
		}

		#Instatiate software instances
		cutadapt = Software('cutadapt', pipeline_config['cutadapt']['path'])
		fastqc = Software('FastQC', pipeline_config['fastqc']['path'])
		bwa_aln = Software('BWA aln', pipeline_config['bwa']['path'] + ' aln')
		bwa_sampe = Software('BWA sampe', pipeline_config['bwa']['path'] + ' sampe')
		samtools_view = Software('samtools view', pipeline_config['samtools']['path'] + ' view')
		samtools_flagstat = Software('samtools flagstat', pipeline_config['samtools']['path'] + ' flagstat')
		samtools_index = Software('samtools index', pipeline_config['samtools']['path'] + ' index')
		samtools_sort = Software('samtools sort', pipeline_config['samtools']['path'] + ' sort')
		novosort = Software('novosort', pipeline_config['novosort']['path'])
		picard_mark_dup = Software('Picard MarkDuplicates', pipeline_config['picard']['path'] + ' MarkDuplicates')
		picard_insert_metrics = Software('Picard CollectInsertSizeMetrics', pipeline_config['picard']['path'] + ' CollectInsertSizeMetrics')
		bedtools_bamtobed = Software('bedtools bamtobed', pipeline_config['bedtools']['path'] + ' bamtobed')
		bedtools_sort = Software('bedtools sort', pipeline_config['bedtools']['path'] + 'sort')
		bedtools_merge = Software('bedtools merge', pipeline_config['bedtools']['path'] + ' merge')
		bedtools_intersect = Software('bedtools intersect', pipeline_config['bedtools']['path'] + ' intersect')
		macs2_callpeak = Software('macs2 callpeak', pipeline_config['macs2']['path'] + ' callpeak')

		if step <= 1:
			for i, read_pair in enumerate(read_pairs):
				read1, read2 = read_pair.split(':')

				#QC: Get raw fastq read counts 
				qc_data['total_raw_reads_counts'].append([
					str(int(self.count_gzipped_lines(read1))/4),
					str(int(self.count_gzipped_lines(read2))/4)
				])

				trimmed_read1_filename = os.path.join(output_dir,
														lib_prefix + '_{}_read1.trimmed.fastq.gz'.format(i))
				trimmed_read2_filename = os.path.join(output_dir,
														lib_prefix + '_{}_read2.trimmed.fastq.gz'.format(i))

				cutadapt.run(
					Parameter('--quality-base=33'),
					Parameter('--minimum-length=5'),
					Parameter('-q',  '30'), # Minimum quality score
					Parameter('--output={}'.format(trimmed_read1_filename)),
					Parameter('--paired-output={}'.format(trimmed_read2_filename)),
					Parameter('-a', forward_adapter if forward_adapter else 'ZZZ'),
					Parameter('-A', reverse_adapter if reverse_adapter else 'ZZZ'),
					Parameter(read1),
					Parameter(read2),
					Redirect(stream=Redirect.STDOUT, dest=os.path.join(logs_dir, 'cutadapt.summary.log'))
				)

				# QC: Get trimmed fastq read counts
				qc_data['trimmed_reads_counts'].append([
					str(int(self.count_gzipped_lines(trimmed_read1_filename))/4),
					str(int(self.count_gzipped_lines(trimmed_read2_filename))/4)
					])

				staging_delete.extend([trimmed_read1_filename, trimmed_read2_filename])
				read_pairs[i] = ':'.join([trimmed_read1_filename, trimmed_read2_filename])

		if step <= 2:
			#Make FastQC Directory
			fastqc_output_dir = os.path.join(output_dir, 'fastqc')
			subprocess.call(['mkdir', '-p', fastqc_output_dir])
			for i, read_pair in enumerate(read_pairs):
				for read in read_pair.split(':'):
					fastqc.run(
						Parameter('--outdir={}'.format(fastqc_output_dir)),
						Parameter(read)
					)

					bwa_aln.run(
						Parameter('-t', pipeline_config['bwa']['threads']),
						Parameter(pipeline_config['bwa']['index-dir']),
						Parameter(read),
						Redirect(stream=Redirect.STDOUT, dest='{}.sai'.format(read))
					)

					staging_delete.append('{}.sai'.format(read))

		if step <= 3:
			for i, read_pair in enumerate(read_pairs):
				read1, read2 = read_pair.split(':')
				bwa_bam_output = os.path.join(output_dir, '{}.{}.bam'.format(lib_prefix, i))

				bwa_sampe.run(
					Parameter('-a', '2000'), # Maximum insert size
					Parameter('-n', '1'),
					Parameter(pipeline_config['bwa']['index-dir']),
					Parameter('{}.sai'.format(read1)),
					Parameter('{}.sai'.format(read2)),
					Parameter(read1),
					Parameter(read2),
					Redirect(stream=Redirect.STDERR, dest=os.path.join(logs_dir, 'bwa_sampe.log')),
					Pipe(
						samtools_view.pipe(
							Parameter('-hSb'),
							Parameter('-o', bwa_bam_output),
							Parameter('-') # Get input from stdin
						)
					)
				)

				bwa_bam_outs.append(bwa_bam_output)

		if step <= 4:
			for i, bwa_bam in enumerate(bwa_bam_outs):
				samtools_flagstat.run(
					Parameter(bwa_bam),
					Redirect(stream=Redirect.STDOUT, dest=bwa_bam + '.flagstat')
				)

				#QC: Get number of mapped reads from this bam
				try:
					with open(bwa_bam + '.flagstat') as flagstats:
						flagstats_contents = flagstats.read()
						target_line = re.search(r'(\d+) \+ \d+ mapped', flagstats_contents)
						if target_line is not None:
							qc_data['num_reads_mapped'].append(str(int(target_line.group(1))/2))
						else:
							qc_data['num_reads_mapped'].append('0')
				except:
					qc_data['num_reads_mapped'].append('Could not open flagstats {}'.format(
						bwa_bam + '.flagstat'
					))

			sortmerged_bam = os.path.join(output_dir, '{}.sortmerged_bam'.format(lib_prefix))
			steric_filter_bam = os.path.join(output_dir, '{}.steric.bam'.format(lib_prefix))
			duprm_bam = os.path.join(output_dir, '{}.duprm.bam'.format(lib_prefix))
			unique_bam = os.path.join(output_dir, '{}.unique.bam'.format(lib_prefix))
			unmappedrm_bam = os.path.join(output_dir, '{}.unmappedrm.bam'.format(lib_prefix))
			chrmrm_bam = os.path.join(output_dir, '{}.chrmrm.bam'.format(lib_prefix))
			# binning read based off template size
			nucleosome_free_reads = os.path.join(output_dir, '{}.nucleosome_free.bam'.format(lib_prefix))
			mononucleosome_reads = os.path.join(output_dir, '{}.mononucleosome.bam'.format(lib_prefix))
			dinucleosome_reads = os.path.join(output_dir, '{}.dinucleosome.bam'.format(lib_prefix))
			trinucleosome_reads = os.path.join(output_dir, '{}.trinucleosome.bam'.format(lib_prefix))
			chrM_bam = os.path.join(output_dir, '{}.chrM.bam'.format(lib_prefix))
			sorted_for_PE_bam = os.path.join(output_dir, '{}.sorted_for_PE'.format(lib_prefix))

			novosort.run(
				Parameter('--threads', pipeline_config['novosort']['threads']),
				Parameter('--tmpcompression', '6'),
				Parameter('--tmpdir', tmp_dir),
				Parameter(*[bam for bam in bwa_bam_outs]),
				Redirect(stream=Redirect.STDOUT, dest=sortmerged_bam),
				Redirect(stream=Redirect.STDERR, dest=os.path.join(logs_dir, 'novosort.log'))
			)

			# This creates a dependency on pysam
			# Removes reads with template length < 38 due to steric hinderence
			samtools_index.run(Parameter(sortmerged_bam))
			sortmerged_bam_alignmentfile = pysam.AlignmentFile(sortmerged_bam, 'rb')
			steric_filter_bam_alignmentfile = pysam.AlignmentFile(steric_filter_bam, 'wb',
																	template=sortmerged_bam_alignmentfile)
			
			num_removed=0
			for read in sortmerged_bam_alignmentfile.fetch():
				if abs(int(read.template_length)) >= STERIC_HINDRANCE_CUTOFF:
					steric_filter_bam_alignmentfile.write(read)
				else:
					num_removed += 1
			qc_data['num_read_removed_steric_hinderence']=str(num_removed)
			
			
			sortmerged_bam_alignmentfile.close()
			steric_filter_bam_alignmentfile.close()

			# Mark and remove MarkDuplicates
			markduplicates_metrics_filepath = os.path.join(logs_dir, 'mark_dup.metrics')
			picard_mark_dup.run(
				Parameter('INPUT={}'.format(steric_filter_bam)),
				Parameter('OUTPUT={}'.format(duprm_bam)),
				Parameter('TMP_DIR={}'.format(tmp_dir)),
				Parameter('METRICS_FILE={}'.format(markduplicates_metrics_filepath)),
				Parameter('REMOVE_DUPLICATES=true'),
				Parameter('VALIDATION_STRINGENCY=LENIENT'),
				Redirect(stream=Redirect.BOTH, dest=os.path.join(logs_dir, 'mark_dup.log'))
			)

			#QC: Get percent MarkDuplicates
			try:
				with open(markduplicates_metrics_filepath) as markdup_metrics:
					for line in markdup_metrics:
						if line[FIRST_CHAR] == '#':
							continue
						record = line.strip().split('\t')
						if len(record) == 9:
							if re.match(r'\d+', record[7]) is not None:
								qc_data['percent_duplicate_reads'] = record[7]
			except:
				qc_data['percent_duplicate_reads'] = 'Could not open MarkDuplicates metrics'

			# Filter down to uniquely mapped reads
			samtools_view.run(
				Parameter('-b'),
				Parameter('-F', '256'),
				Parameter('-q', '10'),
				Parameter('-o', unique_bam),
				Parameter(duprm_bam)
			)

			# gets statistics on uniquely mapped reads
			for i, unique_map in enumerate(unique_bam):
				samtools_flagstat.run(
					Parameter(unique_bam),
					Redirect(stream=Redirect.STDOUT, dest=unique_bam + '.flagstat')
				)

				#QC: Get number of mapped reads from unique bams
			try:
				with open(unique_bam + '.flagstat') as flagstats:
					unique_flagstats_contents = flagstats.read()
					target_line = re.search(r'(\d+) \+ \d+ mapped', unique_flagstats_contents)
					if target_line is not None:
						qc_data['num_unique_reads_mapped'].append(str(int(target_line.group(1))/2))
					else:
						qc_data['num_unique_reads_mapped'].append('0')
			except:
				qc_data['num_unique_reads_mapped'] + '.flagstat'

			# make AlignmentFile object to extract binned reads and chrM reads from the unique bam
			samtools_index.run(Parameter(unique_bam))
			unique_bam_alignmentfile = pysam.AlignmentFile(unique_bam, 'rb')
			# Bins reads into 4 categories depending on template length read is derived from:
			# 50-115 (nucleosome-free), 180-247 (mononucleosome), 315-473 (dinucleosome), 558-615 (trinucleosome)
			nucleosome_free_reads_alignmentfile = pysam.AlignmentFile(nucleosome_free_reads, 'wb',
																	template=unique_bam_alignmentfile)
			mononucleosome_reads_alignmentfile = pysam.AlignmentFile(mononucleosome_reads, 'wb',
																	template=unique_bam_alignmentfile)
			dinucleosome_reads_alignmentfile = pysam.AlignmentFile(dinucleosome_reads, 'wb',
																	template=unique_bam_alignmentfile)
			trinucleosome_reads_alignmentfile = pysam.AlignmentFile(trinucleosome_reads, 'wb',
																	template=unique_bam_alignmentfile)
			
			# Extract chrM into new BAM
			chrM_reads_alignmentfile = pysam.AlignmentFile(chrM_bam, 'wb',
														template=unique_bam_alignmentfile)

			# Binning of nucleosome reads
			for read in unique_bam_alignmentfile.fetch():
				if abs(int(read.template_length)) >= 50 and abs(int(read.template_length)) <= 115:
					nucleosome_free_reads_alignmentfile.write(read)
				elif abs(int(read.template_length)) >= 180 and abs(int(read.template_length)) <= 247:
					mononucleosome_reads_alignmentfile.write(read)
				elif abs(int(read.template_length)) >= 315 and abs(int(read.template_length)) <= 473:
					dinucleosome_reads_alignmentfile.write(read)
				elif abs(int(read.template_length)) >= 558 and abs(int(read.template_length)) <= 615:
					trinucleosome_reads_alignmentfile.write(read)
				else:
					continue;

			#stores chrM reads in separate file
			for read in unique_bam_alignmentfile.fetch():
				if read.reference_name == 'chrM':
					chrM_reads_alignmentfile.write(read)
	
			nucleosome_free_reads_alignmentfile.close()
			mononucleosome_reads_alignmentfile.close()
			dinucleosome_reads_alignmentfile.close()
			trinucleosome_reads_alignmentfile.close()
			chrM_reads_alignmentfile.close()
			
			# gets series of flagstats results for non-main files
			samtools_flagstat.run(
					Parameter(nucleosome_free_reads),
					Redirect(stream=Redirect.STDOUT, dest=nucleosome_free_reads + '.flagstat'))

			samtools_flagstat.run(
					Parameter(mononucleosome_reads),
					Redirect(stream=Redirect.STDOUT, dest=mononucleosome_reads + '.flagstat'))

			samtools_flagstat.run(
					Parameter(dinucleosome_reads),
					Redirect(stream=Redirect.STDOUT, dest=dinucleosome_reads + '.flagstat'))

			samtools_flagstat.run(
					Parameter(trinucleosome_reads),
					Redirect(stream=Redirect.STDOUT, dest=trinucleosome_reads + '.flagstat'))

			
			# gets statistics on chrM mapped reads
			samtools_index.run(Parameter(chrM_bam))
			for i, chrM_map in enumerate(chrM_bam):
				samtools_flagstat.run(
					Parameter(chrM_bam),
					Redirect(stream=Redirect.STDOUT, dest=chrM_bam + '.flagstat')
				)
			try:
				with open(chrM_bam + '.flagstat') as flagstats:
					chrM_flagstats_contents = flagstats.read()
					target_line = re.search(r'(\d+) \+ \d+ mapped', chrM_flagstats_contents)
					if target_line is not None:
						qc_data['num_mtDNA_reads_mapped'].append(str(int(target_line.group(1))/2))
					else:
						qc_data['num_mtDNA_reads_mapped'].append('0')
			except:
				qc_data['num_mtDNA_reads_mapped'] + '.flagstat'



			# Remove unmapped reads
			samtools_view.run(
				Parameter('-b'),
				Parameter('-F', '12'),
				Parameter('-o', unmappedrm_bam),
				Parameter(unique_bam)
			)

			# Create BAM index, then remove chrM
			samtools_index.run(
				Parameter(unmappedrm_bam)
			)

			# Remove chrM
			all_chr = [Parameter('chr{}'.format(chromosome)) for chromosome in map(str, range(1, 23)) + ['X', 'Y']]
			samtools_view.run(
				Parameter('-b'),
				Parameter('-o', chrmrm_bam),
				Parameter(unmappedrm_bam),
				*all_chr
			)

			# Stage delete for temporary files
			staging_delete.extend([
				sortmerged_bam,
				sortmerged_bam + '.bai', # BAM index file
				steric_filter_bam,
				unique_bam,
				duprm_bam,
				unmappedrm_bam,
				unmappedrm_bam + '.bai', # BAM index file
				chrmrm_bam
			])

		if step <= 5:
			# Generate filename for final processed BAM and BED
			processed_bam = os.path.join(output_dir, '{}.processed.bam'.format(lib_prefix))
			unshifted_bed = os.path.join(output_dir, '{}.unshifted_bed'.format(lib_prefix))
			processed_bed = os.path.join(output_dir, '{}.processed.bed'.format(lib_prefix))
			unshifted_bedpe = os.path.join(output_dir, '{}.unshifted_bedpe'.format(lib_prefix))
			processed_bedpe_to_bed = os.path.join(output_dir,'{}.processed_bedpe_to_bed'.format(lib_prefix))
			# staging_delete.append(unshifted_bed)

			# Generate filename for chrM removed BAM
			chrmrm_bam = os.path.join(output_dir, '{}.chrmrm.bam'.format(lib_prefix))

			# Remove blacklisted genomic regions
			bedtools_intersect.run(
				Parameter('-v'),
				Parameter('-abam', chrmrm_bam),
				Parameter('-b', pipeline_config['bedtools']['blacklist-bed']),
				Parameter('-f', '0.5'),
				Redirect(stream=Redirect.STDOUT, dest=processed_bam)
			)

			# QC: Generate insert size metrics PDF
			picard_insert_metrics.run(
				Parameter('INPUT={}'.format(processed_bam)),
				Parameter('OUTPUT={}'.format(os.path.join(logs_dir, lib_prefix + '.insertsize.metrics'))),
				Parameter('HISTOGRAM_FILE={}'.format(os.path.join(logs_dir, lib_prefix + '.insertsize.pdf')))
			)

			# Generate index for processed BAM
			samtools_index.run(
				Parameter(processed_bam)
			)

			# Convert BAM to BED
			bedtools_bamtobed.run(
				Parameter('-i', processed_bam),
				Redirect(stream=Redirect.STDOUT, dest=unshifted_bed)
			)

			# Convert BAM to BEDPE, with specific quality and only properly paired reads, sorted by name
			samtools_view.run(
				Parameter('-uf', '0x2'),
				Parameter('-F', '1548'),
				Parameter('-q', '30'),
				Parameter(processed_bam),
				Pipe(
					samtools_sort.pipe(
						Parameter('-n'),
						Parameter('-'),
						Parameter(sorted_for_PE_bam)
					)
				)
			)

			# convert bam to BEDPE
			bedtools_bamtobed.run(
				Parameter('-i', str(sorted_for_PE_bam)+'.bam'),
				Parameter('-bedpe'),
				Redirect(stream=Redirect.STDOUT, dest=unshifted_bedpe)
			)
			
			unshifted_bedpe_to_bed = open(output_dir+'/'+'{}.unshifted_bedpe_to_bed'.format(lib_prefix), 'w')
			
			with open(unshifted_bedpe) as convertToBed:
				for line in convertToBed:
					chrpos1, start1, end1, chrpos2, start2, end2, name, score, strand1, strand2=line.split('\t')
					bedformat=[chrpos1, start1, end2, name, score, strand1, strand2.rstrip('\n')]
					unshifted_bedpe_to_bed.write('\t'.join(bedformat)+'\n')

					
			staging_delete.append(unshifted_bed)
			staging_delete.append(output_dir+'/'+'{}.unshifted_bedpe_to_bed'.format(lib_prefix))

			# Shifting + strand by 4 and - strand by -5, according to the ATACseq paper

			# This ysed to be bedtools shift, but they are fired
			self.shift_reads(
				input_bed_filepath=unshifted_bed,
				output_bed_filepath=processed_bed,
				log_filepath=os.path.join(logs_dir, 'shift_reads.logs'),
				genome_sizes_filepath=pipeline_config['bedtools']['genome-sizes'],
				minus_strand_shift=MINUS_STRAND_SHIFT,
				plus_strand_shift=PLUS_STRAND_SHIFT
			)

			##TO DO, needs modification for bedpe format
			self.shift_reads_bedpe(
				input_bed_filepath=output_dir+'/'+'{}.unshifted_bedpe_to_bed'.format(lib_prefix),
				output_bed_filepath=processed_bedpe_to_bed,
				log_filepath=os.path.join(logs_dir, 'shift_reads_bedpe_to_bed.logs'),
				genome_sizes_filepath=pipeline_config['bedtools']['genome-sizes'],
				minus_strand_shift=MINUS_STRAND_SHIFT,
				plus_strand_shift=PLUS_STRAND_SHIFT
			)

		# Peak-calling; MACS2
		if step <= 6:
			# for regular peak calling, including narrow, default q-value=0.01
			macs2_callpeak.run(
				Parameter('-t', processed_bed),
				Parameter('-f', 'BED'),
				Parameter('-g', 'hs'),
				Parameter('-n', str(processed_bed) + '_regular_peak_calls'),
				Parameter('--nomodel'),
				Parameter('--extsize', '200'),
				Parameter('--shift', '-100'),
				Parameter('-B', '--SPMR'), # Generates pileup tracks, bedgraph, fragment pileup per million reads
				Parameter('--call-summits'),
				Parameter('--keep-dup', 'all')
			)

			#for broad peak calling, q-value=0.05 per MACS2 suggestion on broad peaks
			macs2_callpeak.run(
				Parameter('-t', processed_bed),
				Parameter('-f', 'BED'),
				Parameter('-g', 'hs'),
				Parameter('-n', str(processed_bed) + '_broad_peak_calls'),
				Parameter('-q', '0.05'),
				Parameter('--nomodel'),
				Parameter('--extsize', '200'),
				Parameter('--shift', '-100'),
				Parameter('--broad'),
				Parameter('--keep-dup', 'all')
			)

			# for regular peak calling, including narrow, default q-value=0.01 for processed bedpe to bed file
			# NOTE: BEDPE for MACS2 is not the same format at BEDPE accepted by NGS/UCSC standards
			macs2_callpeak.run(
				Parameter('-t', processed_bedpe_to_bed),
				Parameter('-f', 'BEDPE'),
				Parameter('-g', 'hs'),
				Parameter('-n', str(processed_bedpe_to_bed) + '_regular_peak_calls'),
				Parameter('--nomodel'),
				Parameter('--extsize', '200'),
				Parameter('--shift', '-100'),
				Parameter('-B', '--SPMR'), # Generates pileup tracks, bedgraph, fragment pileup per million reads
				Parameter('--call-summits'),
				Parameter('--keep-dup', 'all')
			)

			#for broad peak calling, q-value=0.05 per MACS2 suggestion on broad peaks for processed bedpe to bed file
			# NOTE: BEDPE for MACS2 is not the same format at BEDPE accepted by NGS/UCSC standards
			macs2_callpeak.run(
				Parameter('-t', processed_bedpe_to_bed),
				Parameter('-f', 'BEDPE'),
				Parameter('-g', 'hs'),
				Parameter('-n', str(processed_bedpe_to_bed) + '_broad_peak_calls'),
				Parameter('-q', '0.05'),
				Parameter('--nomodel'),
				Parameter('--extsize', '200'),
				Parameter('--shift', '-100'),
				Parameter('--broad'),
				Parameter('--keep-dup', 'all')
			)


		# QC: Output QC data to file
		with open(os.path.join(logs_dir, 'qc_metrics.txt'), 'w') as qc_data_file:
			qc_data_file.write(str(qc_data) + '\n')

		# Delete temporary files
		for delete_file in staging_delete:
			subprocess.call(['rm', '-rf', delete_file])
