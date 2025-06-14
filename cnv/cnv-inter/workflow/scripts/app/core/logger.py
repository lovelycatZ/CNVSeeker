"""
Configure uvicorn logs with loguru for FastAPI
"""
import logging
import http
# import contextlib
from pprint import pformat
from typing import Union

from loguru import logger
from loguru._defaults import LOGURU_FORMAT


class InterceptHandler(logging.Handler):
    """
    Default handler from examples in loguru documentaion.
    See https://loguru.readthedocs.io/en/stable/overview.html#entirely-compatible-with-standard-logging
    """

    def emit(self, record: logging.LogRecord):
        # Get corresponding Loguru level if it exists

        # 处理 uvicorn.access 的内容
        if record.name == 'uvicorn.access':
            record = self.access_format(record)

        try:
            level = logger.level(record.levelname).name
        except ValueError:
            level = record.levelno

        # Find caller from where originated the logged message
        frame, depth = logging.currentframe(), 2
        while frame.f_code.co_filename == logging.__file__:
            frame = frame.f_back
            depth += 1

        logger.opt(depth=depth, exception=record.exc_info, colors=True).log(
            level, record.getMessage()
        )
    
    def access_format(self, record: logging.LogRecord):
        """
        恢复 uvicorn.access message 的默认格式
        默认的 uvicorn 的 access 的 handler 有对 message 做颜色和格式的处理, 直接覆盖 logger 会导致 format 消失

        record 部分内容: (See https://docs.python.org/3/library/logging.html#logging.LogRecord for more details.)
            - record.name: uvicorn.access
            - record.msg:  %s - "%s %s HTTP/%s" %d
            - record.args: ('127.0.0.1:58294', 'GET', '/api/variant/6_162201165_162201165_C_T', '1.1', 200)
            - record.exc_info: None
            - record.exc_text: None
            - record.getMessage(): 127.0.0.1:58294 - "GET /api/variant/6_162201165_162201165_C_T HTTP/1.1" 200
        """
        
        status_code = record.args[-1]
        # 设置描述
        try:
            status_phrase = http.HTTPStatus(status_code).phrase
        except ValueError:
            status_phrase = ""
        status_and_phrase = f"{status_code} {status_phrase}"

        # 设置颜色
        # loguru 设置颜色: See https://loguru.readthedocs.io/en/stable/api/logger.html#loguru._logger.Logger
        if status_and_phrase.startswith("1"):
            status_and_phrase = f'<white>{status_and_phrase}</white>'
        elif status_and_phrase.startswith("2"):
            status_and_phrase = f'<green>{status_and_phrase}</green>'
        elif status_and_phrase.startswith("3"):
            status_and_phrase = f'<yellow>{status_and_phrase}</yellow>'
        elif status_and_phrase.startswith("4"):
            status_and_phrase = f'<light-red>{status_and_phrase}</light-red>'
        elif status_and_phrase.startswith("5"):
            status_and_phrase = f'<red>{status_and_phrase}</red>'

        record.args = list(record.args)
        record.args[-1] = status_and_phrase
        record.args = tuple(record.args)
        record.msg = '%s - "%s %s HTTP/%s" %s'

        return record


class CustomizeLogger:

    @classmethod
    def init_logger(cls, config: Union[dict, str, None] = None):
        """
        Replaces logging handlers with a handler for using the custom handler.
            
        WARNING!
        if you call the init_logging in startup event function, 
        then the first logs before the application start will be in the old format
        >>> app.add_event_handler("startup", init_logging)
        stdout:
        INFO:     Uvicorn running on http://127.0.0.1:8000 (Press CTRL+C to quit)
        INFO:     Started reloader process [11528] using statreload
        INFO:     Started server process [6036]
        INFO:     Waiting for application startup.
        2020-07-25 02:19:21.357 | INFO     | uvicorn.lifespan.on:startup:34 - Application startup complete.
        
        """
        
        # disable handlers for specific uvicorn loggers
        # to redirect their output to the default uvicorn logger
        # works with uvicorn==0.11.6
        loggers = (
            logging.getLogger(name)
            for name in logging.root.manager.loggerDict
            if name.startswith("uvicorn.") # 从 logger 中获取所有 uvicorn 的子 logger 记录器
        )
        for uvicorn_logger in loggers:
            uvicorn_logger.handlers = []    # 清空 logger 的 handlers 处理器

        # change handler for default uvicorn logger
        intercept_handler = InterceptHandler()
        # logging.basicConfig(handlers=[intercept_handler])
        logging.getLogger("uvicorn").handlers = [intercept_handler]
        logging.getLogger("uvicorn.access").handlers = [intercept_handler] # uvicorn.access 的 propagate == False, 不会传递 log record, 因此需要单独处理
        
        # set logs output, level and format
        if isinstance(config, dict):
            logger.configure( **config )

        # # 捕获标准输出流到 logger
        # stream = StreamToLogger()
        # with contextlib.redirect_stdout(stream):
        #     print("Standard output is sent to added handlers.")


class StreamToLogger:

    def __init__(self, level="INFO"):
        self._level = level

    def write(self, buffer):
        for line in buffer.rstrip().splitlines():
            logger.opt(depth=1).log(self._level, line.rstrip())

    def flush(self):
        pass


def format_record(record: dict) -> str:
    """
    Custom format for loguru loggers.
    Uses pformat for log any data like request/response body during debug. 【为了展示请求/相应 body 的内容】
    Works with logging if loguru handler it.
    Example:
    >>> payload = [{"users":[{"name": "Nick", "age": 87, "is_active": True}, {"name": "Alex", "age": 27, "is_active": True}], "count": 2}]
    >>> logger.bind(payload=).debug("users payload")
    >>> [   {   'count': 2,
    >>>         'users': [   {'age': 87, 'is_active': True, 'name': 'Nick'},
    >>>                      {'age': 27, 'is_active': True, 'name': 'Alex'}]}]
    """

    format_string = LOGURU_FORMAT
    if record["extra"].get("payload") is not None:
        record["extra"]["payload"] = pformat(
            record["extra"]["payload"], indent=4, compact=True, width=88
        )
        format_string += "\n<level>{extra[payload]}</level>"

    format_string += "{exception}\n"
    return format_string
