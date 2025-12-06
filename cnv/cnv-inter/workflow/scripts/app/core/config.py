import sys, os
import logging

from typing import Literal
from pydantic import BaseModel, computed_field

from loguru._defaults import LOGURU_FORMAT

from .logger import format_record


def env(key, type_=str, default=None):
    """获取环境变量, 并转换成指定类型"""

    if key not in os.environ:
        return default

    val = os.environ[key]

    if type_ == str:
        return val
    elif type_ == bool:
        if val.lower() in ["1", "true", "yes", "y", "ok", "on"]:
            return True
        if val.lower() in ["0", "false", "no", "n", "nok", "off"]:
            return False
        raise ValueError(
            "Invalid environment variable '%s' (expected a boolean): '%s'" % (key, val)
        )
    elif type_ == int:
        try:
            return int(val)
        except ValueError:
            raise ValueError(
                "Invalid environment variable '%s' (expected an integer): '%s'" % (key, val)
            ) from None
        

class Settings(BaseModel):
    
    PROJECT_NAME: str = 'CNVSeeker'  ####################### 换成你自己的 #####################

    API_STR: str = '/api'
    DOMAIN: str = "localhost"
    ENVIRONMENT: Literal["local", "staging", "production"] = env("_API_ENV", str) if env("_API_ENV") else "local"

    @computed_field  # type: ignore[misc]
    @property
    def app_logger(self) -> dict: 
        logger_config = {
            "handlers": [
                # {"sink": sys.stdout, "level": logging.DEBUG if self.ENVIRONMENT=='local' else logging.INFO, "format": format_record},
                {"sink": sys.__stdout__, "level": logging.DEBUG if self.ENVIRONMENT=='local' else logging.INFO, "format": format_record},
                {
                    "sink": os.path.join(os.path.dirname(__file__), "../logs/app/out.log"),  ####################### 日志位置换成你自己的 #####################
                    "level": logging.DEBUG if self.ENVIRONMENT=='local' else logging.INFO,
                    "format": format_record,
                    "enqueue": True,
                    "backtrace": True,
                    "diagnose": True if self.ENVIRONMENT=='local' else False,
                    "rotation": "1 days",
                    "retention": "1 months"
                },
                {
                    "sink": os.path.join(os.path.dirname(__file__), "../logs/app/err.log"),   ####################### 日志位置换成你自己的 #####################
                    "level": logging.WARNING,
                    "format": format_record,
                    "enqueue": True,
                    "backtrace": True,
                    "diagnose": True if self.ENVIRONMENT=='local' else False,
                    "rotation": "7 days",
                    "retention": "1 months"
                }
            ] 
        }
        return logger_config

    @computed_field  # type: ignore[misc]
    @property
    def server_host(self) -> str:
        # Use HTTPS for anything other than local development
        if self.ENVIRONMENT == "local":
            return f"http://{self.DOMAIN}"
        return f"https://{self.DOMAIN}"



settings = Settings()
