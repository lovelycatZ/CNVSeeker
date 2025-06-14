def get_final_out():
	out_list = []
	if CALL_OUTPUT_FORMAT in ["bam", "cram"]:
		for sample in SAMPLES:
			platform = DF.loc[sample, 'platform']
			if sample in fa_samples + fa_gz_samples:
				out_bam = final_dir + f"{sample}.{platform}.minimap2-asm.{CALL_OUTPUT_FORMAT}"
				out_list.append(out_bam)
			elif sample in fq_samples + fq_gz_samples:
				out_bam = final_dir + f"{sample}.{platform}.minimap2.{CALL_OUTPUT_FORMAT}"
				out_list.append(out_bam)
			else:
				out_bam = final_dir + f"{sample}.{platform}.{CALL_CALL_OUTPUT_FORMAT}"
				out_list.append(out_bam)		
	else:
		if INTEGRATE:
			out_dir = merged_vcf_dir + "output/"
		else:
			out_dir = final_dir

		for sample in SAMPLES:
			out_file = out_dir + f"{sample}.{CALL_OUTPUT_FORMAT}"
			out_list.append(out_file)
	
	return out_list


def get_platform(wildcards):
	sample = wildcards.sample
	platform = str(DF.loc[sample, 'platform'])
	platform = platform.lower()
	return platform


def get_aligner_params(wildcards, aligner):
	params_dict = {
		"minimap2"	: {
			"ont"		: "-x map-ont", 
			"pb-hifi"	: "-x map-hifi",
			"pb-clr"	: "-x map-pb"
		},
		"winnowmap"	: {
			"ont"		: "-x map-ont", 
			"pb-hifi"	: "-x map-pb",
			"pb-clr"	: "-x map-pb-clr"
		},
		"ngmlr"		: {
			"ont"		: "-x ont", 
			"pb-hifi"	: "-x pacbio",
			"pb-clr"	: "-x pacbio"
		},
		"pbmm2"		: {
			"pb-hifi" 	: "--preset HIFI",
			"pb-clr"	: ""
		}
	}
	platform = wildcards.platform
	return params_dict[aligner][platform]


def get_caller_filter_params(tool, option_value_dict):
	caller_filter_options_dict = {
		"cuteSV"		: 	["--min_support", "--min_size", "--max_size"],
		"SVIM"		  	: 	["--minimum_depth", "--min_sv_size", "--max_sv_size"],
		"pbsv"		  	: 	["--call-min-reads-one-sample", "--min-sv-length", ""],
		"SVision_pro" 	: 	["--min_supp", "--min_sv_size", "--max_sv_size"],
		"Sniffles"	  	: 	["", "--minsvlen", ""],
		"NanoVar" 		: 	["--mincov", "--minlen", ""],
		"DeBreak" 		: 	["--min_support", "--min_size", "--max_size"],
		"SVIM_asm"		: 	["", "--min_sv_size", "--max_sv_size"],
		"SVDSS"			:	["--min-cluster-weight", "--min-sv-length", ""]
	}
	params_list = []
	for k, v in zip(caller_filter_options_dict[tool], option_value_dict.values()):
		if k and v:
			params_list.append(f"{k} {v}")
	return " ".join(params_list)
	

def get_caller_input_sample(wildcards, tool):
	tool_aligner_dict = {
		"cuteSV"		: 	"minimap2",
		"SVIM"		  	: 	"minimap2",
		"pbsv"		  	: 	"pbmm2",
		"SVision_pro" 	: 	"minimap2",
		"Sniffles"	  	: 	"minimap2",
		"NanoVar" 		: 	"minimap2",
		"DeBreak" 		: 	"minimap2",
		"SVIM_asm"		: 	"minimap2-asm",
		"SVDSS"			:	"minimap2"
	}
	sample = wildcards.sample
	platform = get_platform(wildcards)
	aligner = tool_aligner_dict[tool]
	if sample in fq_samples + fq_gz_samples + ccs_bam_samples + fa_samples + fa_gz_samples:
		return bam_dir + f"{sample}.{platform}.{aligner}.ready.bam"
	if sample in bam_samples + cram_samples:
		return bam_dir + f"{sample}.{platform}.ready.bam"


def get_caller_output_vcfs(wildcards):
	sample = wildcards.sample
	vcfs = []
	if sample in fa_samples + fa_gz_samples:
		tools = config["used_callers"][DATA_TYPE]["assembly_based"]
	else:
		tools = config["used_callers"][DATA_TYPE]["align_based"]

	for tool in tools:
		vcf = tools_vcf_dir + f"{tool}/output/{sample}.vcf"
		vcfs.append(vcf)
	return vcfs


def get_cuteSV_params(wildcards):
	platform = get_platform(wildcards)
	if platform == 'ont':
		return "--min_read_len 500 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3"
	if platform == 'pb-clr':
		return "--min_read_len 500 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5"
	if platform == 'pb-hifi':
		return "--min_read_len 500 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5"


def get_NanoVar_params(wildcards):
	platform = get_platform(wildcards)
	if platform == 'ont':
		return "-x ont"
	if platform == 'pb-clr':
		return "-x pacbio-clr"
	if platform == 'pb-hifi':
		return "-x pacbio-ccs"


def get_Nanovar_out_vcf(wildcards, file_dir):
	sample = wildcards.sample
	platform = get_platform(wildcards)
	if sample in fq_samples + fq_gz_samples + ccs_bam_samples + fa_samples + fa_gz_samples:
		return file_dir + f"{sample}/{sample}.{platform}.minimap2.ready.nanovar.pass.vcf"
	if sample in bam_samples + cram_samples:
		return file_dir + f"{sample}/{sample}.{platform}.ready.nanovar.pass.vcf"


def get_fastq_sample(wildcards):
	sample = wildcards.sample
	platform = get_platform(wildcards)
	if DF.loc[sample, 'file'].endswith((".fq", ".fastq")):
		return fastq_dir + f"{sample}.{platform}.fq"
	if DF.loc[sample, 'file'].endswith((".fq.gz", ".fastq.gz")) or DF.loc[sample, 'file'].endswith(".ccs.bam"):
		return fastq_dir + f"{sample}.{platform}.fq.gz"


def get_fasta_sample(wildcards):
	sample = wildcards.sample
	platform = get_platform(wildcards)
	if DF.loc[sample, 'file'].endswith((".fa", ".fasta")):
		return assmb_fasta_dir + f"{sample}.{platform}.assmb.fa"
	if DF.loc[sample, 'file'].endswith((".fa.gz", ".fasta.gz")):
		return assmb_fasta_dir + f"{sample}.{platform}.assmb.fa.gz"


rule link_fq:
	input:
		fq = lambda wildcards: DF.loc[DF.index == wildcards.sample, 'file']
	output:
		ln_fq = fastq_dir + "{sample}.{platform}.fq"
	wildcard_constraints:
		sample = "|".join(fq_samples) if fq_samples else "no_wildcards"
	benchmark: benchmark_dir + "link_fq/{sample}.{platform}.log"
	run:
		if not os.path.islink(output[0]):
			os.symlink(input[0], output[0])


rule link_fq_gz:
	input:
		fq = lambda wildcards: DF.loc[DF.index == wildcards.sample, 'file']
	output:
		ln_fq = fastq_dir + "{sample}.{platform}.fq.gz"
	wildcard_constraints:
		sample = "|".join(fq_gz_samples) if fq_gz_samples else "no_wildcards"
	benchmark: benchmark_dir + "link_fq_gz/{sample}.{platform}.log"
	run:
		if not os.path.islink(output[0]):
			os.symlink(input[0], output[0])


rule link_assmb_fa:
	input:
		fa = lambda wildcards: DF.loc[DF.index == wildcards.sample, 'file']
	output:
		ln_fa = assmb_fasta_dir + "{sample}.{platform}.assmb.fa"
	benchmark: benchmark_dir + "link_assmb_fa/{sample}.{platform}.log"
	run:
		if not os.path.islink(output[0]):
			os.symlink(input[0], output[0])


rule link_assmb_fa_gz:
	input:
		fa = lambda wildcards: DF.loc[DF.index == wildcards.sample, 'file']
	output:
		ln_fa = assmb_fasta_dir + "{sample}.{platform}.assmb.fa.gz"
	benchmark: benchmark_dir + "link_assmb_fa_gz/{sample}.{platform}.log"
	run:
		if not os.path.islink(output[0]):
			os.symlink(input[0], output[0])


rule link_bam:
	input:
		bam = lambda wildcards: DF.loc[DF.index == wildcards.sample, 'file'],
		bai = lambda wildcards: DF.loc[DF.index == wildcards.sample, 'file'] + ".bai"
	output:
		ln_bam = bam_dir + "{sample}.{platform}.ready.bam",
		ln_bai = bam_dir + "{sample}.{platform}.ready.bam.bai"
	benchmark: benchmark_dir + "link_bam/{sample}.{platform}.log"
	run:
		if not os.path.islink(output[0]):
			os.symlink(input[0], output[0])
		if not os.path.islink(output[1]):
			os.symlink(input[1], output[1])


rule link_cram:
	input:
		cram = lambda wildcards: DF.loc[DF.index == wildcards.sample, 'file'],
		crai = lambda wildcards: DF.loc[DF.index == wildcards.sample, 'file'] + ".crai"
	output:
		ln_cram = bam_dir + "{sample}.{platform}.ready.cram",
		ln_crai = bam_dir + "{sample}.{platform}.ready.cram.crai"
	benchmark: benchmark_dir + "link_cram/{sample}.{platform}.log"
	run:
		if not os.path.islink(output[0]):
			os.symlink(input[0], output[0])
		if not os.path.islink(output[1]):
			os.symlink(input[1], output[1])


rule cram2bam:
	input:
		cram = bam_dir + "{sample}.{platform}.ready.cram",
		crai = bam_dir + "{sample}.{platform}.ready.cram.crai"
	output:
		bam = bam_dir + "{sample}.{platform}.ready.bam",
		bai = bam_dir + "{sample}.{platform}.ready.bam.bai"
	params:
		ref = REFERENCE
	conda: TGS_env_1
	benchmark: benchmark_dir + "cram2bam/{sample}.{platform}.log"
	log: logs_dir + "cram2bam/{sample}.{platform}.log"
	threads: THREADS
	shell:
		"(samtools view -Sb -T {params.ref} -@ {threads} {input.cram} > {output.bam}; "
		"samtools index -@ {threads} {output.bai} {output.bai}; "
		") > {log} 2>&1"


rule link_ccs_bam:
	input:
		ccs_bam = lambda wildcards: DF.loc[DF.index == wildcards.sample, 'file'],
		pbi = lambda wildcards: DF.loc[DF.index == wildcards.sample, 'file'] + ".pbi"
	output:
		ln_ccs_bam = ccs_bam_dir + "{sample}.{platform}.ccs.bam",
		ln_pbi = ccs_bam_dir + "{sample}.{platform}.ccs.bam.pbi"
	benchmark: benchmark_dir + "link_ccs_bam/{sample}.{platform}.log"
	run:
		if not os.path.islink(output[0]):
			os.symlink(input[0], output[0])
		if not os.path.islink(output[1]):
			os.symlink(input[1], output[1])


rule ccs2fastq:
	input:
		ccs_bam = ccs_bam_dir + "{sample}.{platform}.ccs.bam"
	output:
		fastq = fastq_dir + "{sample}.{platform}.fq.gz"
	params:
		out_prefix = fastq_dir + "{sample}.{platform}",
		fastq = fastq_dir + "{sample}.{platform}.fastq.gz"
	conda: TGS_env_1
	benchmark: benchmark_dir + "ccs2fastq/{sample}.{platform}.log"
	log: logs_dir + "ccs2fastq/{sample}.{platform}.log"
	threads: THREADS
	wildcard_constraints:
		sample = "|".join(ccs_bam_samples) if ccs_bam_samples else "no_wildcards"
	shell:
		"(bam2fastq -j {threads} -o {params.out_prefix} {input.ccs_bam}; "
		"mv {params.fastq} {output.fastq}; "
		") > {log} 2>&1"
