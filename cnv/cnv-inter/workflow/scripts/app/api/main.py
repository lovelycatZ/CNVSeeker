from fastapi import APIRouter

from app.api.routes import cnv
from app.api.routes import file


api_router = APIRouter()

api_router.include_router(cnv.router, prefix='/cnv')
api_router.include_router(file.router, prefix='/file')
