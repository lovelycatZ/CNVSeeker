from fastapi import APIRouter, BackgroundTasks, status

from app.services.file import interpret_file
from app.schemas.acmg import FileReqBody


router = APIRouter()

@router.post("", response_model=dict, status_code=status.HTTP_202_ACCEPTED)
def get_file_report(file_in: FileReqBody, background_tasks: BackgroundTasks, genome_build:str='hg38'):
    
    file_path = file_in.path
    outdir, outname = file_in.out_dir, file_in.out_name
    callback = file_in.callback

    background_tasks.add_task(interpret_file, file_path, outdir, outname, genome_build=genome_build, callback=callback)

    return {"message": "Task submitted successfully!"}
