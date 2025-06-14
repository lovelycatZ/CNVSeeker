from typing import Optional
from pydantic import BaseModel


class ServiceError(BaseModel):

    msg: str
    traceback: Optional[str] = None
