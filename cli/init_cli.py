import argparse

import argcomplete

from cli.cmd.cnv_cmd import CNVCallCmd, CNVInterCmd, CNVCallInterCmd
from cli.cmd.conf_helper import ConfHelperCmd

from cli.init_data import CNV_CALL, CNV_INTER, CNV_CALL_INTER
from cli.init_data import CONF_HELPER
from cli.args.cnv_args import cnv_call_args, cnv_inter_args, cnv_call_inter_args
from cli.args.conf_helper import config_helper_args


subcommands = [
    CNVCallCmd,
    CNVInterCmd,
    CNVCallInterCmd,
    ConfHelperCmd,
]

def init_parser(subcommands, des=None) -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    # 构建子命令
    for command in subcommands:  ## type: Type[Command]
        if command.name == CNV_CALL:
            cmd_instance = command(cnv_call_args)
        elif command.name == CNV_INTER:
            cmd_instance = command(cnv_inter_args)
        elif command.name == CNV_CALL_INTER:
            cmd_instance = command(cnv_call_inter_args)
        elif command.name == CONF_HELPER:
            cmd_instance = command(config_helper_args)
        else:
            ...
        # 添加子命令
        subparser = subparsers.add_parser(command.name) 
        # 设置默认处理方法
        subparser.set_defaults(handle=cmd_instance.handle) # handle 参数的名字自行定义

        cmd_instance.add_group_arguments(subparser) # 添加分组参数 (无法添加分组时，会自行调用普通添加函数)

    argcomplete.autocomplete(parser)

    return parser.parse_args()
