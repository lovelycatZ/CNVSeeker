## Automatically update home directory of cnvseeker

## once should provide resource dir and conda env dir, just run bash configure.sh!
resource_dir=/gpfs/hpc/home/lijc/xiangxud/project/resources
conda_env_dir=/gpfs/hpc/home/lijc/xiangxud/software/miniconda/envs
##


resource_dir=${resource_dir%/}
conda_env_dir=${conda_env_dir%/}


cnvseeker_home_name=BASE_DIR
cnv_call_inter_home_name=cnv_call_inter_home
cnv_call_home_name=cnv_call_home
cnv_inter_home_name=cnv_inter_home

cnvseeker_home=$(pwd)
cnv_call_inter_home=$cnvseeker_home/cnv/cnv-call-inter
cnv_call_home=$cnvseeker_home/cnv/cnv-call
cnv_inter_home=$cnvseeker_home/cnv/cnv-inter


# update cnvseeker config
sed -i 's#^'$cnvseeker_home_name': .*#'$cnvseeker_home_name': '$cnvseeker_home'#' $cnvseeker_home/config/config.yml

# update cnv call and interpretation config
sed -i 's#^'$cnv_call_inter_home_name': .*#'$cnv_call_inter_home_name': '$cnv_call_inter_home'#' $cnv_call_inter_home/config/config.yaml
sed -i 's#^'$cnv_call_home_name': .*#'$cnv_call_home_name': '$cnv_call_home'#' $cnv_call_inter_home/config/config.yaml
sed -i 's#^'$cnv_inter_home_name': .*#'$cnv_inter_home_name': '$cnv_inter_home'#' $cnv_call_inter_home/config/config.yaml
sed -i 's#^resource_dir: .*#resource_dir: '$resource_dir'#' $cnv_call_inter_home/config/config.yaml
sed -i 's#^\(.*conda_env_dir\).*#\1: '$conda_env_dir'#' $cnv_call_inter_home/config/config.yaml

# update cnv call config
sed -i 's#^'$cnv_call_home_name': .*#'$cnv_call_home_name': '$cnv_call_home'#' $cnv_call_home/config/config.yaml
sed -i 's#^resource_dir: .*#resource_dir: '$resource_dir'#' $cnv_call_home/config/config.yaml
sed -i 's#^\(.*conda_env_dir\).*#\1: '$conda_env_dir'#' $cnv_call_home/config/config.yaml

# update cnv interpretation config
sed -i 's#^'$cnv_inter_home_name': .*#'$cnv_inter_home_name': '$cnv_inter_home'#' $cnv_inter_home/config/config.yaml
sed -i 's#^resource_dir: .*#resource_dir: '$resource_dir'#' $cnv_inter_home/config/config.yaml

# update cnv interpretation resource
annotation_dir=${resource_dir}/annotation_datasets
sed -i 's#^resource_dir =.*#resource_dir = "'$annotation_dir'"#' $cnv_inter_home/workflow/scripts/interpretation/resources.py

# update autopvs1 config
autopvs1_dir=${resource_dir}/annotation_datasets/autopvs1_data
example_config=${cnv_inter_home}/workflow/scripts/interpretation/autopvs1/config.example.ini
new_config=${cnv_inter_home}/workflow/scripts/interpretation/autopvs1/config.ini
sed "s|data/|${autopvs1_dir}/|g" "$example_config" > "$new_config"


# Tab 补全功能激活
activate-global-python-argcomplete --user
