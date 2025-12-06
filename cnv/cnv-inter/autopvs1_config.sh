resource_dir=/gpfs/hpc/home/lijc/xiangxud/project/resources
conda_env_path=/gpfs/hpc/home/lijc/xiangxud/software/miniconda3/envs
ref_hg38=${resource_dir}/ref_genome/hg38/hg38.fa
ref_hg19=${resource_dir}/ref_genome/hg19/hg19.fa

vep_r=112
cache_r=104
vep_cache_37=${resource_dir}/vep_caches/homo_sapiens_refseq_vep_${cache_r}_GRCh37.tar.gz
vep_cache_38=${resource_dir}/vep_caches/homo_sapiens_refseq_vep_${cache_r}_GRCh38.tar.gz


conda_env_name=interpret_env
python_version=3.11
vep_path=${conda_env_path}/${conda_env_name}/share/ensembl-vep-${vep_r}.0-0
autopvs1_config=${conda_env_path}/${conda_env_name}/lib/python${python_version}/site-packages/autopvs1/config.ini


## vep config
# FTP='ftp://ftp.ensembl.org/pub/'
# r=112
# indexed vep cache
cd $vep_path
# wget -b $FTP/release-${cache_r}/variation/indexed_vep_cache/homo_sapiens_refseq_vep_${cache_r}_GRCh38.tar.gz & \
cp $vep_cache_37 ./
tar xzf homo_sapiens_refseq_vep_${cache_r}_GRCh37.tar.gz
# wget -b $FTP/release-${cache_r}/variation/indexed_vep_cache/homo_sapiens_refseq_vep_${cache_r}_GRCh37.tar.gz & \
cp $vep_cache_38 ./
tar xzf homo_sapiens_refseq_vep_${cache_r}_GRCh38.tar.gz

# fasta
source activate $conda_env_name
cd ${vep_path}/homo_sapiens_refseq/${cache_r}_GRCh37/
ln -s $ref_hg19 hg19.fa
samtools faidx hg19.fa
cd ${vep_path}/homo_sapiens_refseq/${cache_r}_GRCh38/
ln -s $ref_hg38 hg38.fa
samtools faidx hg38.fa


## autpvs1 config
sed -i "2s#.*#vep_cache = $vep_path#" $autopvs1_config
sed -i "8s#^.*#genome = $ref_hg19#" $autopvs1_config
sed -i "17s#^.*#genome = $ref_hg38#" $autopvs1_config