#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK
from cli.init_cli import init_parser, subcommands


def main():
    args = init_parser(subcommands=subcommands)

    # 处理子命令对应的函数
    # 需要和 set_defaults 一致
    if hasattr(args, 'handle'):
        args.handle(args)


if __name__ == '__main__':

    main()