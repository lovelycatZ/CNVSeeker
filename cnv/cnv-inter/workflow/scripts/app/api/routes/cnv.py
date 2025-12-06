from fastapi import APIRouter, BackgroundTasks, status

from app.services.cnv import interpret_cnv_record, interpret_pedigree_cnv_record, interpret_cnv_records
from app.schemas.acmg import VariantsReqBody, VariantReqBody

router = APIRouter()

@router.get("/")
async def root():
    return {"message": "Hello CNVSeeker"}


@router.get("/{cnv_str}", response_model = dict)
def get_cnv_info(cnv_str: str, genome_build: str = "hg38"):

    res = interpret_cnv_record(cnv_str, genome_build)
    
    return res


@router.post("/cnv_pedigree_str", response_model = dict)
def get_pedigree_cnv_info(var_req_body: VariantReqBody, genome_build: str = "hg38"):

    pedigree_cnv_str = var_req_body.cnv_str
    relation_list = var_req_body.relation_list
    res = interpret_pedigree_cnv_record(pedigree_cnv_str, genome_build, relation_list)
    
    return res


@router.post("/multiple_records", response_model=dict, status_code=status.HTTP_202_ACCEPTED)
def get_cnv_records(vars_req_body: VariantsReqBody, background_tasks: BackgroundTasks, genome_build:str='hg38'):
    
    cnv_list = vars_req_body.cnv_list
    outdir, outname = vars_req_body.out_dir, vars_req_body.out_name
    callback = vars_req_body.callback

    background_tasks.add_task(interpret_cnv_records, cnv_list, outdir, outname, genome_build=genome_build, callback=callback)

    return {"message": "Task submitted successfully!"}
