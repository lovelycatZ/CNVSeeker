from interpretation.utils import DataBase


resource_dir = "/gpfs/hpc/home/lijc/xiangxud/project/resources/annotation_datasets"

chain_file_T2T = resource_dir + "/chain_files/T2TTohg38.chain.gz"
file_dir_hg19 = resource_dir + "/hg19/"
file_dir_hg38 = resource_dir + "/hg38/"

# hg19 soure files
gene_info_hg19 = file_dir_hg19 + "genes.info.hg19.sort.bed.gz"
cytoband_hg19 = file_dir_hg19 + "cytoBand_hg19.sort.bed.gz"
ClinGen_regions_hg19 = file_dir_hg19 + "ClinGen_region_curation_list_GRCh37.sort.bed.gz"
ClinVar_hg19 = file_dir_hg19 + "clinvar.hg19.sort.bed.gz"
ClinVar_null_variant_hg19 = file_dir_hg19 + "clinvar_null_variant.hg19.sort.bed.gz"
gnomAD_SV_hg19 = file_dir_hg19 + "gnomad.v2.1_sv.sites.hg19.sort.bed.gz"
DGV_GS_CNV_hg19 = file_dir_hg19 + "DGV_GS_hg19.sort.bed.gz"
DDD_CNV_hg19 = file_dir_hg19 + "DECIPHER_population_cnv.hg19.sort.bed.gz"
ClinVar_CNV_hg19 = file_dir_hg19 + "clinvar_CNV.hg19.sort.bed.gz"
#dbvar_common_hg19 = file_dir_hg19 + "dbvar.nr_common_CNVs.GRCh37.sort.bed.gz"
#dbvar_pathogenic_hg19 = file_dir_hg19 + "dbvar.nr_pathogenic_CNVs.GRCh37.sort.bed.gz"

# hg38 soure files
gene_info_hg38 = file_dir_hg38 + "genes.info.hg38.sort.bed.gz"
cytoband_hg38 = file_dir_hg38 + "cytoBand_hg38.sort.bed.gz"
ClinGen_regions_hg38 = file_dir_hg38 + "ClinGen_region_curation_list_GRCh38.sort.bed.gz"
ClinVar_hg38 = file_dir_hg38 + "clinvar.hg38.sort.bed.gz"
ClinVar_null_variant_hg38 = file_dir_hg38 + "clinvar_null_variant.hg38.sort.bed.gz"
gnomAD_SV_hg38 = file_dir_hg38 + "gnomad.v4.1_sv.sites.hg38.sort.bed.gz"
DGV_GS_CNV_hg38 = file_dir_hg38 + "DGV_GS_hg38.sort.bed.gz"
DDD_CNV_hg38 = file_dir_hg38 + "DECIPHER_population_cnv.hg38.sort.bed.gz"
ClinVar_CNV_hg38 = file_dir_hg38 + "clinvar_CNV.hg38.sort.bed.gz"
#dbvar_common_hg38 = file_dir_hg38 + "dbvar.nr_common_CNVs.GRCh38.sort.bed.gz"
#dbvar_pathogenic_hg38 = file_dir_hg38 + "dbvar.nr_pathogenic_CNVs.GRCh38.sort.bed.gz"


# hg19 databases
DB_gene_info_hg19 = DataBase(gene_info_hg19)
DB_cytoband_hg19 = DataBase(cytoband_hg19)
DB_ClinGen_regions_hg19 = DataBase(ClinGen_regions_hg19)
DB_ClinVar_hg19 = DataBase(ClinVar_hg19)
DB_ClinVar_null_variant_hg19 = DataBase(ClinVar_null_variant_hg19)
DB_gnomAD_SV_hg19 = DataBase(gnomAD_SV_hg19)
DB_DGV_GS_CNV_hg19 = DataBase(DGV_GS_CNV_hg19)
DB_DDD_CNV_hg19 = DataBase(DDD_CNV_hg19)
DB_ClinVar_CNV_hg19 = DataBase(ClinVar_CNV_hg19)
#DB_dbvar_common_hg19 = DataBase(dbvar_common_hg19)
#DB_dbvar_pathogenic_hg19 = DataBase(dbvar_pathogenic_hg19)

# hg38 databases
DB_gene_info_hg38 = DataBase(gene_info_hg38)
DB_cytoband_hg38 = DataBase(cytoband_hg38)
DB_ClinGen_regions_hg38 = DataBase(ClinGen_regions_hg38)
DB_ClinVar_hg38 = DataBase(ClinVar_hg38)
DB_ClinVar_null_variant_hg38 = DataBase(ClinVar_null_variant_hg38)
DB_gnomAD_SV_hg38 = DataBase(gnomAD_SV_hg38)
DB_DGV_GS_CNV_hg38 = DataBase(DGV_GS_CNV_hg38)
DB_DDD_CNV_hg38 = DataBase(DDD_CNV_hg38)
DB_ClinVar_CNV_hg38 = DataBase(ClinVar_CNV_hg38)
#DB_dbvar_common_hg38 = DataBase(dbvar_common_hg38)
#DB_dbvar_pathogenic_hg38 = DataBase(dbvar_pathogenic_hg38)
