import logging

from fastapi import FastAPI
import uvicorn

from app.core.config import settings
from app.core.logger import CustomizeLogger
from app.api.main import api_router


app = FastAPI(
    title=settings.PROJECT_NAME,
    openapi_url=f"{settings.API_STR}/openapi.json"
)

# initial logger
CustomizeLogger.init_logger(settings.app_logger)


# Set all CORS enabled origins
# TODO


print(settings)

app.include_router(api_router, prefix=settings.API_STR)
