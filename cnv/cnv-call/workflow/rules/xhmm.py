## 获取覆盖深度
rule get_depth:
	input:
		bam = bam_dir + "{sample}.ready.bam"
	output:
		TOOLS_DIR_DICT["xhmm"]["info"] + "depth/{sample}.sample_interval_summary"
	params:
		outdir = TOOLS_DIR_DICT["xhmm"]["info"] + "depth/{sample}",
		bed = EXOME_INTERVALS,
		ref = REFERENCE
	conda: NGS_env_1
	benchmark: benchmark_dir + "xhmm/get_depth/{sample}.log"
	log: logs_dir + "xhmm/get_depth/{sample}.log"
	threads: THREADS
	shell:
		"(gatk3 -T DepthOfCoverage "
		"-I {input.bam} -L {params.bed} -R {params.ref} "
		"-dt BY_SAMPLE -dcov 5000 -l INFO -omitBaseOutput -omitLocusTable "
		"--minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 "
		"--includeRefNSites "
		"--countType COUNT_FRAGMENTS "
		"-o {params.outdir}; "
		") > {log} 2>&1"


## 合并深度结果
rule merge_depth:
	input:
		expand(TOOLS_DIR_DICT["xhmm"]["info"] + "depth/{sample}.sample_interval_summary", sample = SAMPLES)
	output:
		TOOLS_DIR_DICT["xhmm"]["info"] + "data.rd.txt"
	params:
		depth_file = TOOLS_DIR_DICT["xhmm"]["info"] + "read_depth_file_list"
	conda: NGS_env_1
	benchmark: benchmark_dir + "xhmm/merge_depth.log"
	log: logs_dir + "xhmm/merge_depth.log"
	threads: THREADS
	shell:
		"(ls {input} > {params.depth_file}; "		 
		"xhmm --mergeGATKdepths --GATKdepthsList {params.depth_file} -o {output}; "
		") > {log} 2>&1"


## 获得target区域GC比例
rule get_GC:
	output:
		locus_gc = TOOLS_DIR_DICT["xhmm"]["info"] + "data.locus.gc.txt",
		extre_gc = TOOLS_DIR_DICT["xhmm"]["info"] + "extreme.gc.target.txt"
	params:
		bed = EXOME_INTERVALS,
		ref = REFERENCE
	conda: NGS_env_1
	benchmark: benchmark_dir + "xhmm/get_GC.log"
	log: logs_dir + "xhmm/get_GC.log"
	threads: THREADS
	shell:
		"(gatk3 -T GCContentByInterval -L {params.bed} -R {params.ref} -o {output.locus_gc}; "
		"cat {output.locus_gc} | awk '{{if ($2 < 0.1 || $2 > 0.9) print $1}}' > {output.extre_gc}; "
		") > {log} 2>&1"


## 计算重复碱基并过滤
# rule calculate_repeat_base:
# 	output:
# 		TOOLS_DIR_DICT["xhmm"]["info"] + "low_complexity_targets.txt"
# 	params:
# 		reg = TOOLS_DIR_DICT["xhmm"]["info"] + "targets.reg",
# 		locdb = TOOLS_DIR_DICT["xhmm"]["info"] + "targets.LOCDB",
# 		loc_load = TOOLS_DIR_DICT["xhmm"]["info"] + "targets.LOCDB.loc-load",
# 		locus_complexity = TOOLS_DIR_DICT["xhmm"]["info"] + "locus_complexity.txt",
# 		pseq = config["tools"]["xhmm"]["bin"]["--pseq"],
# 		seqdb = config["tools"]["xhmm"]["params"]["--seqdb"],
# 		bed = EXOME_INTERVALS,
# 	log: logs_dir + "xhmm/calculate_repeat_base.log"
# 	threads: THREADS
# 	shell:
# 		"(echo -e \"#CHR\tBP1\tBP2\tID\" > {params.reg}; "
# 		"awk -v OFS='\t' '{{print $0, NR}}' {params.bed} >> {params.reg}; "
# 		"{params.pseq} . loc-load --locdb {params.locdb} --file {params.reg} --out {params.loc_load} --group targets --noweb; "
# 		"{params.pseq} . loc-stats --locdb {params.locdb} --group targets --seqdb {params.seqdb} | "
# 		"awk '{{if (NR > 1) print $_}}' | sort -k1 -g | awk '{{print $10}}' | paste {params.reg} - > {params.locus_complexity}; "
# 		"cat {params.locus_complexity} | awk '{{if ($5 > 0.25) print $1 \":\" $2 \"-\" $3}}' > {output}; "
# 		") > {log} 2>&1"


## 过滤
rule filter:
	input:
		rd = TOOLS_DIR_DICT["xhmm"]["info"] + "data.rd.txt",
		extreme_gc_target = TOOLS_DIR_DICT["xhmm"]["info"] + "extreme.gc.target.txt",
		# low_complexity_targets = TOOLS_DIR_DICT["xhmm"]["info"] + "low_complexity_targets.txt"
	output:
		centered_rd = TOOLS_DIR_DICT["xhmm"]["info"] + "data.filtered_centered.RD.txt",
		out_targets = TOOLS_DIR_DICT["xhmm"]["info"] + "data.filtered_centered.RD.txt.filtered_targets.txt",
		out_samples = TOOLS_DIR_DICT["xhmm"]["info"] + "data.filtered_centered.RD.txt.filtered_samples.txt"
	conda: NGS_env_1
	benchmark: benchmark_dir + "xhmm/filter.log"
	log: logs_dir + "xhmm/filter.log"
	threads: THREADS
	shell:
		"(xhmm --matrix "
		"-r {input.rd} --centerData --centerType target "
		"-o {output.centered_rd} "
		"--outputExcludedTargets {output.out_targets} "
		"--outputExcludedSamples {output.out_samples} "
		"--excludeTargets {input.extreme_gc_target} "
		# "--excludeTargets {input.low_complexity_targets} "
		"--minTargetSize 10 --maxTargetSize 10000 "
		"--minMeanTargetRD 10 --maxMeanTargetRD 500 "
		"--minMeanSampleRD 25 --maxMeanSampleRD 200 "
		"--maxSdSampleRD 150; "
		") > {log} 2>&1"


## PCA聚类分析
rule PCA:
	input:
		TOOLS_DIR_DICT["xhmm"]["info"] + "data.filtered_centered.RD.txt"
	output:
		TOOLS_DIR_DICT["xhmm"]["info"] + "data.RD.PCA.PC_LOADINGS.txt",
		TOOLS_DIR_DICT["xhmm"]["info"] + "data.RD.PCA.PC_SD.txt",
		TOOLS_DIR_DICT["xhmm"]["info"] + "data.RD.PCA.PC.txt"
	params:
		pca = TOOLS_DIR_DICT["xhmm"]["info"] + "data.RD.PCA"
	conda: NGS_env_1
	benchmark: benchmark_dir + "xhmm/PCA.log"
	log: logs_dir + "xhmm/PCA.log"
	threads: THREADS
	shell:
		"(xhmm --PCA -r {input} --PCAfiles {params.pca}; "
		") > {log} 2>&1"


## 降噪
rule normalize:
	input:
		centered_rd = TOOLS_DIR_DICT["xhmm"]["info"] + "data.filtered_centered.RD.txt",
		PCA_loadings = TOOLS_DIR_DICT["xhmm"]["info"] + "data.RD.PCA.PC_LOADINGS.txt",
		PCA_sd = TOOLS_DIR_DICT["xhmm"]["info"] + "data.RD.PCA.PC_SD.txt",
		PCA_pc = TOOLS_DIR_DICT["xhmm"]["info"] + "data.RD.PCA.PC.txt"
	output:
		TOOLS_DIR_DICT["xhmm"]["info"] + "data.PCA_normalized.txt"
	params:
		pca = TOOLS_DIR_DICT["xhmm"]["info"] + "data.RD.PCA"
	conda: NGS_env_1
	benchmark: benchmark_dir + "xhmm/normalize.log"
	log: logs_dir + "xhmm/normalize.log"
	threads: THREADS
	shell:
		"(xhmm --normalize -r {input.centered_rd} --PCAfiles {params.pca} --normalizeOutput {output} "
		"--PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7; "
		") > {log} 2>&1"


## Z检验
rule Z_test:
	input:
		pca = TOOLS_DIR_DICT["xhmm"]["info"] + "data.PCA_normalized.txt"
	output:
		zscore = TOOLS_DIR_DICT["xhmm"]["info"] + "data.PCA_normalized.filtered.sample_zscores.RD.txt",
		out_targets = TOOLS_DIR_DICT["xhmm"]["info"] + "data.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt",
		out_samples = TOOLS_DIR_DICT["xhmm"]["info"] + "data.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt"
	conda: NGS_env_1
	benchmark: benchmark_dir + "xhmm/Z_test.log"
	log: logs_dir + "xhmm/Z_test.log"
	threads: THREADS
	shell:
		"(xhmm --matrix "
		"-r {input.pca} --centerData --centerType sample --zScoreData "
		"-o {output.zscore} "
		"--outputExcludedTargets {output.out_targets} "
		"--outputExcludedSamples {output.out_samples} "
		"--maxSdTargetRD 30; "
		") > {log} 2>&1"


## 原始深度过滤
rule filter_RD:
	input:
		rd = TOOLS_DIR_DICT["xhmm"]["info"] + "data.rd.txt",
		pca = TOOLS_DIR_DICT["xhmm"]["info"] + "data.PCA_normalized.txt",
		out_targets1 = TOOLS_DIR_DICT["xhmm"]["info"] + "data.filtered_centered.RD.txt.filtered_targets.txt",
		out_samples1 = TOOLS_DIR_DICT["xhmm"]["info"] + "data.filtered_centered.RD.txt.filtered_samples.txt",
		out_targets2 = TOOLS_DIR_DICT["xhmm"]["info"] + "data.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt",
		out_samples2 = TOOLS_DIR_DICT["xhmm"]["info"] + "data.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt"
	output:
		TOOLS_DIR_DICT["xhmm"]["info"] + "data.same_filtered.RD.txt"
	conda: NGS_env_1
	benchmark: benchmark_dir + "xhmm/filter_RD.log"
	log: logs_dir + "xhmm/filter_RD.log"
	threads: THREADS
	shell:
		"(xhmm --matrix "
		"-r {input.rd} "
		"--excludeTargets {input.out_targets1} "
		"--excludeTargets {input.out_targets2} "
		"--excludeSamples {input.out_samples1} "
		"--excludeSamples {input.out_samples2} "
		"-o {output}; "
		") > {log} 2>&1"


## call CNV
rule xhmm_CNV_detect:
	input:
		zscore = TOOLS_DIR_DICT["xhmm"]["info"] + "data.PCA_normalized.filtered.sample_zscores.RD.txt",
		same_rd = TOOLS_DIR_DICT["xhmm"]["info"] + "data.same_filtered.RD.txt"
	output:
		xcnv = TOOLS_DIR_DICT["xhmm"]["info"] + "data.xcnv",
		aux_xcnv = TOOLS_DIR_DICT["xhmm"]["info"] + "data.aux_xcnv"
	params:
		options = config["tools"]["xhmm"]["params"]["--options"]
	conda: NGS_env_1
	benchmark: benchmark_dir + "xhmm/xhmm_CNV_detect.log"
	log: logs_dir + "xhmm/xhmm_CNV_detect.log"
	threads: THREADS
	shell:
		"(xhmm --discover "
		"-p {params.options} "
		"-r {input.zscore} -R {input.same_rd} "
		"-c {output.xcnv} -a {output.aux_xcnv}; "
		") > {log} 2>&1"

		
## 导出为VCF格式
rule export_vcf:
	input:
		ref = REFERENCE,
		zscore = TOOLS_DIR_DICT["xhmm"]["info"] + "data.PCA_normalized.filtered.sample_zscores.RD.txt",
		same_rd = TOOLS_DIR_DICT["xhmm"]["info"] + "data.same_filtered.RD.txt",
		xcnv = TOOLS_DIR_DICT["xhmm"]["info"] + "data.xcnv"
	output:
		vcf = TOOLS_DIR_DICT["xhmm"]["info"] + "data.xhmm.vcf"
	params:
		options = config["tools"]["xhmm"]["params"]["--options"]
	conda: NGS_env_1
	benchmark: benchmark_dir + "xhmm/export_vcf.log"
	log: logs_dir + "xhmm/export_vcf.log"
	threads: THREADS
	shell:
		"(xhmm --genotype "
		"-p {params.options} "
		"-r {input.zscore} -R {input.same_rd} "
		"-g {input.xcnv} -F {input.ref} "
		"-v {output.vcf}; "
		") > {log} 2>&1"


## 转换格式
rule xhmm_postprocess:
	input:
		vcf = TOOLS_DIR_DICT["xhmm"]["info"] + "data.xhmm.vcf"
	output:
		vcf = TOOLS_DIR_DICT["xhmm"]["output"] + "{sample}_xhmm.vcf"
	benchmark: benchmark_dir + "xhmm/xhmm_postprocess/{sample}.log"
	run:
		xhmm_postprocess(input[0], wildcards.sample, output[0])
		# add_var_id(output[0], "xhmm")
