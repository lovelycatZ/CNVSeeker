# CNVSeeker: an analytical pipeline for empowering comprehensive analysis of clinical disease-associated copy number variations from multiple sequencing strategies

`CNVSeeker` is a comprehensive ensemble one-stop pipeline that provided a series of capacities such as data    pre-process, detection, annotation and interpretation of CNVs from raw sequencing data to variant interpretation report. CNVSeeker consists of two main modules: cnv-call (CNV calling) and cnv-inter (CNV interpretation). cnv-call is designed to detect CNVs from next-generation sequencing (NGS) and third-generation sequencing (TGS) while cnv-inter focus on pathogenicity interpretation of germline CNVs. Various type of input formats are permitted such as FASTQ (.gz), BAM, CRAM, SAM, VCF, and BED.

![image](https://github.com/lovelycatZ/CNVSeeker/blob/main/framework.svg)

## Quick Start Guideline

### Step 1: Clone from github
```
git clone https://github.com/lovelycatZ/CNVSeeker.git
```

### Step 2: Download resource data
The resource datasets used for implementation of CNVSeeker have been deposited at Zenodo (https://doi.org/10.5281/zenodo.15660867).
```
wget https://doi.org/10.5281/zenodo.15660867
tar -zxvf /path/to/cnvseeker_resource.tar.gz
```

### Step 3: Configure the necessary properties
You should first provide **resource directory ($resource_dir)** and **conda enviroment directory ($conda_env_dir)** in shell script `configure.sh`, and then just run `bash configure.sh`. This step will automatically update home directory of varseeker and its submodules.
```
cd /path/to/CNVSeeker/cnv/
bash configure.sh
```

### Step 4: Create conda environment

Since too many tools or packages are used in our pipline, we choose to manage them using conda.

**Create conda environment for TGS module:**
```
cd /path/to/CNVSeeker/cnv/cnv-call/workflow/envs/
conda env create -f conda.TGS.1.yml
conda env create -f conda.TGS.2.yml
```

**Create conda environment for NGS module:**
```
cd /path/to/CNVSeeker/cnv/cnv-call/workflow/envs/
conda env create -f conda.NGS.1.yml
conda env create -f conda.NGS.2.yml
conda env create -f ECOLE.yml
conda env create -f gatk4.yml
```

**Create conda environment for snakemake:**
```
cd /path/to/CNVSeeker/cnv/cnv-call/workflow/envs/
conda env create -f snakemake.yml
```

**Create conda environment for CNV interpretation module:**
```
cd /path/to/CNVSeeker/cnv/cnv-inter/workflow/envs/
conda env create -f conda.interpret.yml
```

## Quick run
```
python cnvseeker.py cnv-call-inter \
  --sample-list=sample.list \
  --result-dir=result_dir \
  --data-type=NGS-HD_WGS \
  --build=hg38 \
  --chrom-style=l \
  --smk '-j=5,-c=5,--use-conda'
```

## Full usage

```
conda activate snakemake

cd /path/to/CNVSeeker
python cnvseeker.py -h

there are three subcommand available in CNVSeeekr:

1. for cnv calling 
python cnvseeker.py cnv-call [options]

2. for cnv interpretation
python cnvseeker.py cnv-inter [options]

3. for cnv calling - interpretation
python cnvseeker.py cnv-call-inter [options]
```

whole pipline example:
```
python cnveeker.py cnv-call-inter -h

usage: cnvseeker.py cnv-call-inter [-h] --sample-list SAMPLE_LIST --result-dir RESULT_DIR --data-type
                                   {TGS,NGS-WES,NGS-LD_WGS,NGS-HD_WGS} -b {hg19,hg38,T2T} --chrom-style {l,s}
                                   [--exome-intervals EXOME_INTERVALS] [--call-output-format {cram,bam,vcf,bed}]
                                   [--inter-output-format {txt,tsv,csv}] [-ref REFERENCE] [-t THREADS] [-smk SNAKEMAKE]
                                   [--fastp FASTP] [--bwa BWA] [--sambamba SAMBAMBA] [--survivor SURVIVOR] [--cnvpytor CNVPYTOR]
                                   [--controlfreec CONTROLFREEC] [--ecole ECOLE] [--cn_mops CN_MOPS] [--delly DELLY] [--lumpy LUMPY]
                                   [--xhmm XHMM] [--gatk4 GATK4] [--TGS-caller-filter-options TGS_CALLER_FILTER_OPTIONS]
                                   [--filtlong FILTLONG] [--svision-pro SVISION_PRO] [--sniffles SNIFFLES] [--svdss SVDSS]

options:
  -h, --help            show this help message and exit

base:
  --sample-list SAMPLE_LIST
                        sample list file. Please refer to the notes in sample.list.* for details [mandatory]
  --result-dir RESULT_DIR
                        the output directory of CNVSeeker [mandatory]
  --data-type {TGS,NGS-WES,NGS-LD_WGS,NGS-HD_WGS}
                        sequencing type of NGS data. NGS high depth WGS (NGS-HD_WGS), NGS low depth WGS (NGS-LD_WGS), NGS WES (NGS-
                        WES), or TGS [mandatory]
  -b {hg19,hg38,T2T}, --build {hg19,hg38,T2T}
                        the reference genome version [hg19, hg38, T2T] [mandatory]
  --chrom-style {l,s}   the style of chromosome IDs. 'l' = 'long style' (eg. 'chr1', 'chrX'); 's' = 'short style' (eg. '1', 'X').
                        [mandatory]
  --exome-intervals EXOME_INTERVALS
                        exome capture intervals in bed format for WES data. It must be in bed format and no header (1st column:
                        chromosome number, 2nd column: start site-1, 3rd column: end site) [mandatory if WES]
  --call-output-format {cram,bam,vcf,bed}
                        output file format for CNV calling [default: vcf]
  --inter-output-format {txt,tsv,csv}
                        output file format for CNV interpreting [default: tsv]
  -ref REFERENCE, --reference REFERENCE
                        the reference genome in fasta format [default:
                        /gpfs/hpc/home/lijc/xiangxud/project/resources/call_datasets/ref_genome/hg38/hg38.fa]
  -t THREADS, --threads THREADS
                        number of threads for each rule [default: 8]

snakemake:
  -smk SNAKEMAKE, --snakemake SNAKEMAKE
                        naive arguments provided to snakemake

tools:
  --fastp FASTP         naive arguments of fastp. Please refer to help page of fastp and provided in a string like this: -z=4,--
                        reads_to_process=0,-A
  --bwa BWA
  --sambamba SAMBAMBA
  --survivor SURVIVOR
  --cnvpytor CNVPYTOR
  --controlfreec CONTROLFREEC
  --ecole ECOLE
  --cn_mops CN_MOPS
  --delly DELLY
  --lumpy LUMPY
  --xhmm XHMM
  --gatk4 GATK4
  --TGS-caller-filter-options TGS_CALLER_FILTER_OPTIONS
                        min_support: minimum number of supporting read for a SV to be reported; min_size: minimum SV length to be
                        reported; max_size: maximum SV length to be reported;
  --filtlong FILTLONG
  --svision-pro SVISION_PRO
  --sniffles SNIFFLES
  --svdss SVDSS
```

**All parameters consists of three parts: base, snakemake and tools：**

### base section
You should provide the base information of pipline, such as the final result directory, 
the data type of input data, output format of final results, sample list file which 
contains a list of case or control samples.

**for CNV calling - interpretation pipline required parameters:**
* --sample-list
* --result-dir
* --exome-bed (must provide if data type is NGS-WES)
* --data-type
* --build
* --chrom-style

**for CNV calling pipline required parameters:**
* --sample-list
* --result-dir
* --exome-bed (must provide if data type is NGS-WES)
* --data-type
* --build
* --chrom-style

**for CNV interpretation pipline required parameters:**
* --inter-input-format
* --result-dir
* --build


### snakemake section
You can provide snakemake native parameters.

**for local use:**

`python varseeker.py cnv-inter [base option] -smk '-j=5,-c=5,--use-conda' [tools option]`


**for cluster use:**
- slurm job scheduling system
  
`python cnvseeker.py cnv-inter [base option] -smk '-j=5,--use-conda,--cluster "sbatch -p NODE_NAME --cpus-per-task=5"' [tools option]`

- qsub job scheduling system
  
`python cnvseeker.py cnv-inter [base option] -smk '-j=5,--use-conda,--cluster "qsub -q NODE -l nodes=1:ppn=5"' [tools option]`


### tools section
You can provide native parameters of which callers or tools used in pipline. For example,
if you would like to change parameters of *fastp*, please refer to the help page of *fastp*,
and do like this:

`python cnvseeker.py cnv-call [base option] [snakemake option] --fastp '-z=4,--reads_to_process=0,-A'`

## Contact
If you are trouble in using CNVSeeker, please feel free to contact with Xudong Xiang (lovelyxxyz@163.com) or report a bug on github.

## License
CNVSeeker is free for non-commercial use by academic and non-profit/not-for-profit institutions.
