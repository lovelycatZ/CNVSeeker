from pydantic import BaseModel
from typing import Optional, List



class VariantReqBody(BaseModel):
    
    cnv_str: str
    relation_list: List[dict]

    
class VariantsReqBody(BaseModel):
    
    cnv_list: List[str]
    out_dir: Optional[str] = None
    out_name: Optional[str] = None
    callback: str


class FileReqBody(BaseModel):
    
    path: str
    out_dir: Optional[str] = None
    out_name: Optional[str] = None
    callback: str
