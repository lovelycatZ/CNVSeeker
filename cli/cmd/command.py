import os
from copy import deepcopy

from snakemake.utils import update_config

from cli.init_data import config_loader
from cli.utils import parse_cmd_arguments, reverse_parse_cmd_arguments, convert_config


class Command(object):
    """commamd 基类"""

    name = 'CMD'

    def __init__(self, arguments) -> None:

        self.arguments = arguments
        # print('========= DEBUG Command init ======== \n', arguments)

    def add_arguments(self, parser):
        # 判断传递的参数是否基本合法
        # if type(self.arguments) is not list:
        if not isinstance(self.arguments, list):
            print(f'EEROR: Need to provide Arguments list. Current is : {type(self.arguments)}')
            return
        for arg in self.arguments:
            arg.add_to_parser(parser)
        # self.subcommand_add_arguments(parser)
            
    def add_group_arguments(self, parser):
        # 判断传递的参数是否基本合法
        # if type(self.arguments) is list:
        if isinstance(self.arguments, list):
            # print(f'WARNING: Need to provide Arguments dict. Current is : {type(self.arguments)}')
            # print(f'WARNING: Trying add arguments by "add_arguments" function')
            self.add_arguments(parser)
            return
        # elif type(self.arguments) is not dict:
        elif not isinstance(self.arguments, dict):
            # print(f'EEROR: Need to provide Arguments dict. Current is : {type(self.arguments)}')
            return

        for group_name, args in self.arguments.items():
            group = parser.add_argument_group(group_name)
            for arg in args:
                arg.add_to_group(group)

    def append_arguments(self, parser, arg):
        arg.add_to_parser(parser)

    def handle(self, cur_cmd, args, args_keychain_map, _config, final_config_path, work_dir) -> None:
        print(f"===DEBUG===>> handle ----command args: {args}")

        config = deepcopy(_config)

        self.cur_cmd = cur_cmd

        # 处理 base 的参数
        base_args_map = args_keychain_map['base'] # args_keychain_map 为 argparse 中的参数（变量名）和配置文件中的参数位置的映射
        self._handle_base(vars(args), base_args_map, config=config)

        # 处理 tools 的参数
        tools_args_map = args_keychain_map['tools']
        self._handle_tools(vars(args), tools_args_map, config=config)

        convert_config(config)
        config_loader.file_dump(config, final_config_path)
        
        # 处理 snakemake（最后处理）
        self._handle_snakemake(args=args.snakemake, config=config, config_path=final_config_path, work_dir=work_dir)

    def _handle_snakemake(self, args, config, config_path, work_dir):

        print(f"===DEBUG===>> _handle_snakemake args: {args}")
        # 转换成字典（类似 tools 中的参数字典）
        smk_args_dict = parse_cmd_arguments(args)  # snakemake 命令行参数，不能以 - 开头（配置文件中可以）

        smk_params_dict_default = config['snakemake']['params'] # snakemake 默认参数
        if isinstance(smk_params_dict_default, dict):
            for smk_p, smk_pval in smk_params_dict_default.items():
                if not isinstance(smk_pval, dict):
                    continue
                arg_name = str(smk_pval['arg_name']).strip() if 'arg_name' in smk_pval.keys() else None
                arg_short_name = str(smk_pval['arg_short_name']).strip() if 'arg_short_name' in smk_pval.keys() else None
                if arg_name in smk_args_dict.keys() or arg_short_name in smk_args_dict.keys():
                    continue
                if arg_name: # 存在 arg_name 没有设置的情况
                    # smk_pval 不一定都是 str，可能是 int 等其它类型（配置文件中读取）
                    smk_args_dict[arg_name] = smk_pval['value'].strip() if isinstance(smk_pval['value'], str) else smk_pval['value']
                elif arg_short_name:
                    smk_args_dict[arg_short_name] = smk_pval['value'].strip() if isinstance(smk_pval['value'], str) else smk_pval['value']

        # 配置文件
        if '--configfile' not in smk_args_dict.keys():  # TODO 配置文件还有一个 --configfiles 没有考虑，后续需要优化
            smk_args_dict['--configfile'] = config_path

        # snakemake 主目录
        if '-d' not in smk_args_dict.keys() and '--directory' not in smk_args_dict.keys():
            smk_args_dict['--directory'] = work_dir

        smk_cmd = reverse_parse_cmd_arguments(smk_args_dict)
        if not smk_cmd.startswith('snakemake '):
            smk_cmd = 'snakemake ' + smk_cmd

        # 处理环境变量
        # environs = config['envs']['environ']
        # env_cmd = ''
        # if hasattr(environs, "__iter__"):
        #     for env in environs:
        #         env_cmd += env.strip() + ' && '
        # smk_cmd = env_cmd + smk_cmd
        print('============ DEBUG smk cmd ============\n cmd: ', smk_cmd)
        os.system(smk_cmd)

    def _handle_base(self, args, args_map, config):
        base_args_dict = {}
        for arg, key_chain in args_map.items():
            if args[arg] is not None:
                val = base_args_dict
                for key in key_chain.split('.'):
                    if val.get(key) is None:
                        val[key] = {}
                    val = val[key]
                val['value'] =args[arg]
        # print('======== DEBUG base dict =======\n', base_args_dict)
        update_config(config, base_args_dict)

    def _handle_tools(self, args, args_map, config):

        tools_args_dict = {}
        for arg, key_chain in args_map.items():
            if args[arg] is not None:
                val = tools_args_dict
                for key in key_chain.split('.'):
                    if val.get(key) is None:
                        val[key] = {}
                    val = val[key]
                params_dict = parse_cmd_arguments(args[arg])
                # print('============= DEBUG prams dict ===========\n', params_dict)
                val['params'] = params_dict
        # print('======== DEBUG tools dict =======\n', tools_args_dict)
        update_config(config, tools_args_dict)
