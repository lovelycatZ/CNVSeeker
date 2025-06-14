import os, traceback
import requests

from loguru import logger
from app.schemas.msg import ServiceError

from interpretation.interpret import interpret


def interpret_file(file_path: str, outdir: str, outname: str, genome_build:str='hg38', callback:str=None, is_backgroud: bool = True):
    try:
        print(os.getcwd())
        outdir = outdir if outdir else os.getcwd() + "/web_output"
        if outname:
            outname = outname
        else:
            if ".vcf" in file_path:
                outname = os.path.basename(file_path).split(".vcf")[0] + "_report.tsv"
            else:
                outname = os.path.basename(file_path).split(".bed")[0] + "_report.tsv"  

        if not os.path.isdir(outdir):
            os.makedirs(outdir)
            logger.debug(f"Create directory: <{outdir}>.")
        
        outfile = os.path.join(outdir, outname)

        if os.path.isfile(outfile):
            logger.warning(f"The file <{outfile}> already exists!")
        else:
            logger.info(f"Output file is <{outfile}>.")

        cores = 20
        interpret(file_path, cores, genome_build, outfile)
    
        url = callback + "?code=1"
        response = requests.get(url)
        if response.status_code != 200:
            logger.warning(f"response from {url}: {response.status_code}")

        return outfile
    
    except Exception as e:
        res = ServiceError(
            msg = str(e),
            traceback = traceback.format_exc()
        )
        logger.exception(f"{e}\n")

        if is_backgroud:
            url = callback + "?code=0"
            requests.get(url=url)
        else:
            return res
