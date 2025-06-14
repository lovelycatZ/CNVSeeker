import os

scatter_count = config["tools"]["gatk4"]["params"]["--scatter_count"]
SCATTER_INDEXS = [i for i in range(1, scatter_count + 1)]



if NGS_DATA_TYPE != "WES":
	rule preprocess_intervals:
		input:
			ref = REFERENCE,
		output:
			preprocessed_interval = TOOLS_DIR_DICT["gatk4"]["info"] + "preprocessed.interval_list"
		params:
			window_size = config["tools"]["gatk4"]["params"]["--window_size"]
		conda: gatk4
		threads: THREADS
		benchmark: benchmark_dir + "gatk4/preprocess_intervals.log"
		log: logs_dir + "gatk4/preprocess_intervals.log"
		shell: 
			"(gatk PreprocessIntervals "
			"-R {input.ref} "
			"--padding 0 "
			"--bin-length {params.window_size} "
			"-imr OVERLAPPING_ONLY "
			"-O {output.preprocessed_interval}; "
			") > {log} 2>&1"
else:
	rule preprocess_intervals:
		input:
			ref = REFERENCE,
			interval = EXOME_INTERVALS
		output:
			preprocessed_interval = TOOLS_DIR_DICT["gatk4"]["info"] + "preprocessed.interval_list"
		params:
			padding_size = config["tools"]["gatk4"]["params"]["--padding_size"]
		conda: gatk4
		threads: THREADS
		benchmark: benchmark_dir + "gatk4/preprocess_intervals.log"
		log: logs_dir + "gatk4/preprocess_intervals.log"
		shell: 
			"(gatk PreprocessIntervals "
			"-R {input.ref} "
			"-L {input.interval} "
			"--padding {params.padding_size} "
			"--bin-length 0 "
			"-imr OVERLAPPING_ONLY "
			"-O {output.preprocessed_interval}; "
			") > {log} 2>&1"	

rule remove_intervals:
	input:
		preprocessed_interval = TOOLS_DIR_DICT["gatk4"]["info"] + "preprocessed.interval_list"
	output:
		preprocessed_new_interval = TOOLS_DIR_DICT["gatk4"]["info"] + "preprocessed.new.interval_list"
	params:
		chrom_list = CHROM_LIST
	benchmark: benchmark_dir + "gatk4/remove_intervals.log"
	run:
		chrom_list = params["chrom_list"]
		input_file = input["preprocessed_interval"]
		output_file = output["preprocessed_new_interval"]
		with open(input_file, "r") as f, open(output_file, "w") as out:
			for line in f:
				if line.startswith("@"):
					out.write(line)
					continue
				chr = line.strip().split("\t")[0]
				if chr not in chrom_list:
					continue
				out.write(line)
		

rule annotate_intervals:
	input:
		preprocessed_new_interval = TOOLS_DIR_DICT["gatk4"]["info"] + "preprocessed.new.interval_list"
	output:
		annotated_interval = TOOLS_DIR_DICT["gatk4"]["info"] + "preprocessed.annotated.tsv",
	params:
		ref = REFERENCE
	conda: gatk4
	threads: THREADS
	benchmark: benchmark_dir + "gatk4/annotate_intervals.log"
	log: logs_dir + "gatk4/annotate_intervals.log"
	shell:
		"(gatk AnnotateIntervals "
		"-L {input.preprocessed_new_interval} "
		"-R {params.ref} "
		"-imr OVERLAPPING_ONLY "
		"-O {output.annotated_interval}; "
		") > {log} 2>&1"
		

rule collect_read_counts:
	input:
		bam = bam_dir + "{sample}.ready.bam",	
		preprocessed_new_interval = TOOLS_DIR_DICT["gatk4"]["info"] + "preprocessed.new.interval_list"
	output:
		read_count = temp(TOOLS_DIR_DICT["gatk4"]["info"] + "read_counts/{sample}.hdf5"),
	params:
		ref = REFERENCE
	conda: gatk4
	threads: THREADS
	benchmark: benchmark_dir + "gatk4/collect_read_counts/{sample}.log"
	log: logs_dir + "gatk4/collect_read_counts/{sample}.log"
	priority: 1
	shell:
		"(gatk CollectReadCounts "
		"-L {input.preprocessed_new_interval} "
		"-R {params.ref} "
		"-imr OVERLAPPING_ONLY "
		"-I {input.bam} "
		"--format HDF5 "
		"-O {output.read_count} "
		"--verbosity DEBUG; "
		") > {log} 2>&1"


rule filter_intervals:
	input:
		read_counts = expand(TOOLS_DIR_DICT["gatk4"]["info"] + "read_counts/{sample}.hdf5", sample = SAMPLES),
		annotated_interval = TOOLS_DIR_DICT["gatk4"]["info"] + "preprocessed.annotated.tsv"
	output:
		gc_content_interval = TOOLS_DIR_DICT["gatk4"]["info"] + "gc.filtered.interval_list"
	params:
		read_count_str = create_read_count_str(file_dir = TOOLS_DIR_DICT["gatk4"]["info"] + "read_counts/", suffix = "hdf5"),
		preprocessed_new_interval = TOOLS_DIR_DICT["gatk4"]["info"] + "preprocessed.new.interval_list"
	conda: gatk4
	threads: THREADS
	benchmark: benchmark_dir + "gatk4/filter_intervals.log"
	log: logs_dir + "gatk4/filter_intervals.log"
	shell:
		"(gatk FilterIntervals "
		"-L {params.preprocessed_new_interval} "
		"--annotated-intervals {input.annotated_interval} "
		"{params.read_count_str} "
		"-imr OVERLAPPING_ONLY "
		"-O {output.gc_content_interval}; "
		") > {log} 2>&1"


rule determine_ploidy:
	input:
		read_counts = expand(TOOLS_DIR_DICT["gatk4"]["info"] + "read_counts/{sample}.hdf5", sample = SAMPLES),
		gc_content_interval = TOOLS_DIR_DICT["gatk4"]["info"] + "gc.filtered.interval_list"
	output:
		ploidy_calls = directory(temp(TOOLS_DIR_DICT["gatk4"]["info"] + "ploidy/ploidy-calls")),
		ploidy_model = directory(temp(TOOLS_DIR_DICT["gatk4"]["info"] + "ploidy/ploidy-model")),
		theano_determine_ploidy = directory(temp(TOOLS_DIR_DICT["gatk4"]["info"] + "theano_flags_determine_ploidy"))
	params:
		outdir = TOOLS_DIR_DICT["gatk4"]["info"] + "ploidy",
		read_count_str = create_read_count_str(file_dir = TOOLS_DIR_DICT["gatk4"]["info"] + "read_counts/", suffix = "hdf5"),
		ploidy_priors = config["tools"]["gatk4"]["params"]["--ploidy_priors_long_chr"] if CHROM_STYLE == "l" else config["tools"]["gatk4"]["params"]["--ploidy_priors_short_chr"],
		prefix = "ploidy",
	conda: gatk4
	threads: THREADS
	benchmark: benchmark_dir + "gatk4/determine_germline_contig_ploidy.log"
	log: logs_dir + "gatk4/determine_germline_contig_ploidy.log"
	shell:
		"(export PYTENSOR_FLAGS='base_compiledir={output.theano_determine_ploidy}'; "
		"gatk DetermineGermlineContigPloidy "
		"-L {input.gc_content_interval} "
		"--interval-merging-rule OVERLAPPING_ONLY "
		"{params.read_count_str} "
		"--contig-ploidy-priors {params.ploidy_priors} "
		"--output {params.outdir} "
		"--output-prefix {params.prefix} "
		"--verbosity DEBUG; "
		") > {log} 2>&1"


rule scatter:
	input:
		gc_content_interval = TOOLS_DIR_DICT["gatk4"]["info"] + "gc.filtered.interval_list"
	output:
		scatter_dir = directory(TOOLS_DIR_DICT["gatk4"]["info"] + "scatter")
	params:
		scatter_count = scatter_count,
	conda: gatk4
	threads: THREADS
	benchmark: benchmark_dir + "gatk4/scatter.log"
	log: logs_dir + "gatk4/scatter.log"
	shell:
		"(mkdir {output.scatter_dir}; "
		"gatk IntervalListTools "
		"--INPUT {input.gc_content_interval} "
		"--SUBDIVISION_MODE INTERVAL_COUNT "
		"--SCATTER_COUNT {params.scatter_count} "
		"--OUTPUT {output.scatter_dir}; "
		") > {log} 2>&1"
 

rule GATK_gCNV_call:
	input:
		scatter_dir = TOOLS_DIR_DICT["gatk4"]["info"] + "scatter",
		read_counts = expand(TOOLS_DIR_DICT["gatk4"]["info"] + "read_counts/{sample}.hdf5", sample = SAMPLES),
		ploidy_calls = TOOLS_DIR_DICT["gatk4"]["info"] + "ploidy/ploidy-calls",
	output:
		calls = directory(temp(TOOLS_DIR_DICT["gatk4"]["info"] + "cohort_all/cohort_{index}-calls")),
		model = directory(temp(TOOLS_DIR_DICT["gatk4"]["info"] + "cohort_all/cohort_{index}-model")),
		tracking = directory(temp(TOOLS_DIR_DICT["gatk4"]["info"] + "cohort_all/cohort_{index}-tracking")),
		theano_germline_caller = directory(temp(TOOLS_DIR_DICT["gatk4"]["info"] + "theano_flags_germline_caller_{index}"))
	params:
		outdir = TOOLS_DIR_DICT["gatk4"]["info"] + "cohort_all",
		read_count_str = create_read_count_str(file_dir = TOOLS_DIR_DICT["gatk4"]["info"] + "read_counts/", suffix = "hdf5"),
		annotated = TOOLS_DIR_DICT["gatk4"]["info"] + "preprocessed.annotated.tsv",
		interval = lambda wildcards: f"{TOOLS_DIR_DICT['gatk4']['info']}scatter/temp_{(4 - len(wildcards.index)) * '0' + wildcards.index}_of_{scatter_count}/scattered.interval_list",
		prefix = "cohort_{index}"
	conda: gatk4
	threads: THREADS
	benchmark: benchmark_dir + "gatk4/germline_CNVCaller/{index}.log"
	log: logs_dir + "gatk4/germline_CNVCaller/{index}.log"
	shell:
		"(export PYTENSOR_FLAGS='base_compiledir={output.theano_germline_caller}'; "
		"gatk GermlineCNVCaller "
		"--run-mode COHORT "
		"-L {params.interval} "
		"{params.read_count_str} "
		"--contig-ploidy-calls {input.ploidy_calls} "
		"--annotated-intervals {params.annotated} "
		"--interval-merging-rule OVERLAPPING_ONLY "
		"--output {params.outdir} "
		"--verbosity DEBUG "
		"--output-prefix {params.prefix} "
		"--max-advi-iter-first-epoch 2500 "
		"--max-advi-iter-subsequent-epochs 100 "
		"--max-training-epochs 25 "
		"--min-training-epochs 7; "
		# "--class-coherence-length 20000.0 "
		# "--cnv-coherence-length 20000.0 "
		# "--interval-psi-scale 1.0E-6 "
		# "--log-mean-bias-standard-deviation 0.01 "
		# "--sample-psi-scale 1.0E-6; "
		") > {log} 2>&1"


rule postprocess_calls:
	input:
		calls = expand(TOOLS_DIR_DICT["gatk4"]["info"] + "cohort_all/cohort_{index}-calls", index = SCATTER_INDEXS),
		model = expand(TOOLS_DIR_DICT["gatk4"]["info"] + "cohort_all/cohort_{index}-model", index = SCATTER_INDEXS),
		ploidy_calls = TOOLS_DIR_DICT["gatk4"]["info"] + "ploidy/ploidy-calls"
	output:
		vcf = temp(TOOLS_DIR_DICT["gatk4"]["info"] + "output/sample_{sample_index}.cohort.intervals.vcf"),
		vcf_idx = temp(TOOLS_DIR_DICT["gatk4"]["info"] + "output/sample_{sample_index}.cohort.intervals.vcf.idx"),
		segment_vcf = temp(TOOLS_DIR_DICT["gatk4"]["info"] + "output/sample_{sample_index}.cohort.segments.vcf"),
		segment_vcf_idx = temp(TOOLS_DIR_DICT["gatk4"]["info"] + "output/sample_{sample_index}.cohort.segments.vcf.idx"),
		ratio = temp(TOOLS_DIR_DICT["gatk4"]["info"] + "output/sample_{sample_index}.ratio.txt"),
		theano_postprocess_calls = directory(temp(TOOLS_DIR_DICT["gatk4"]["info"] + "theano_flags_postprocess_calls_{sample_index}"))
	params:
		calls_shard_str = create_calls_shard_str(TOOLS_DIR_DICT["gatk4"]["info"] + "cohort_all/", SCATTER_INDEXS),
		model_shard_str = create_model_shard_str(TOOLS_DIR_DICT["gatk4"]["info"] + "cohort_all/", SCATTER_INDEXS),
		allosomal_contig_X = "chrX" if CHROM_STYLE == "l" else "X",
		allosomal_contig_Y = "chrY" if CHROM_STYLE == "l" else "Y"
	conda: gatk4
	threads: THREADS
	benchmark: benchmark_dir + "gatk4/postprocess_germlineCNVCalls/sample_{sample_index}.log"
	log: logs_dir + "gatk4/postprocess_germlineCNVCalls/sample_{sample_index}.log"
	shell:
		"(export PYTENSOR_FLAGS='base_compiledir={output.theano_postprocess_calls}'; "
		"gatk PostprocessGermlineCNVCalls "
		"{params.model_shard_str} {params.calls_shard_str} "
		"--contig-ploidy-calls {input.ploidy_calls} "
		"--allosomal-contig {params.allosomal_contig_X} "
		"--allosomal-contig {params.allosomal_contig_Y} "
		"--sample-index {wildcards.sample_index} "
		"--output-genotyped-intervals {output.vcf} "
		"--output-genotyped-segments {output.segment_vcf} "
		"--output-denoised-copy-ratios {output.ratio}; "
		") > {log} 2>&1"


SAMPLE_INDEXS = [i for i in range(len(SAMPLES))]
rule index2sample:
	input:
		segment_vcfs = expand(TOOLS_DIR_DICT["gatk4"]["info"] + "output/sample_{sample_index}.cohort.segments.vcf", sample_index = SAMPLE_INDEXS),
		ploidy_calls = TOOLS_DIR_DICT["gatk4"]["info"] + "ploidy/ploidy-calls"
	output:
		segment_vcf = TOOLS_DIR_DICT["gatk4"]["info"] + "output/{sample}.cohort.segments.vcf"
	params:
		sample_name = "{sample}"
	threads: THREADS
	benchmark: benchmark_dir + "gatk4/index2sample/{sample}.log"
	run:
		for segment_vcf in input["segment_vcfs"]:
			index = segment_vcf.split("sample_")[-1].split(".")[0]
			sample_file = input["ploidy_calls"] + f"/SAMPLE_{index}/sample_name.txt"
			with open(sample_file, "r") as f:
				sample_name = f.readline().rstrip()
			if sample_name == params["sample_name"]:
				out_segment_vcf = output["segment_vcf"]
				os.system(f"cp {segment_vcf} {out_segment_vcf}")
				break


rule gatk4_postprocess:
	input:
		vcf = TOOLS_DIR_DICT["gatk4"]["info"] + "output/{sample}.cohort.segments.vcf"
	output:
		vcf = protected(TOOLS_DIR_DICT["gatk4"]["output"] + "{sample}_gatk4.vcf")
	benchmark: benchmark_dir + "gatk4/gatk4_postprocess/{sample}.log"
	run:
		gatk4_postprocess(input[0], wildcards.sample, output[0])
		# add_var_id(output[0], "GATK_gCNV")
