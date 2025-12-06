import os
from copy import deepcopy

from cli.utils import convert_config
from cli.cmd.command import Command
from cli.init_data import cmd_config, config_loader, BASE_DIR, CONF_HELPER
from cli.init_data import CNV_CALL, CNV_INTER, CNV_CALL_INTER


class ConfHelperCmd(Command):
    """ CONFIG_HELPER command """

    name = CONF_HELPER

    def __init__(self, arguments) -> None:
        super().__init__(arguments)

        self.name = CONF_HELPER
        self.work_dir = os.path.join(BASE_DIR, self.name)


    def append_arguments(self, parser, arg):
        super().append_arguments(parser, arg)

    def handle(self, args) -> None:
        # print(f"===DEBUG===>> config helper args: {args}")
        if args.create:
            config_user = os.path.join(BASE_DIR, 'config', f"user.{args.create}_config.yml")
            if os.path.isfile(config_user):
                print(f"\nERROR: Configuration file of command <{args.create}> already exists!\n")
                return
            
            self._create_config(args.create)
            print("\nCreate successfully!\n")

        elif args.reset:
            config_user = os.path.join(BASE_DIR, 'config', f"user.{args.reset}_config.yml")
            if not os.path.isfile(config_user):
                print(f"\nERROR: Configuration file of command <{args.reset}> does not exist!\n")
                return
            res = input(f"\nAre you sure you want to reset the configuration file of command <{args.reset}>? (Y/N) ")
            if res.strip().lower() in ['y', 'yes']:
                os.remove(config_user)
                self._create_config(args.reset)
                print("\nReset successfully!\n")


    def _get_default_config(self, cmd_name):
        """根据命令名称获取对应命令的默认配置文件"""

        if cmd_name == CNV_CALL:
            _default_config = cmd_config.cnv_call_config
        elif cmd_name == CNV_INTER:
            _default_config = cmd_config.cnv_inter_config
        elif cmd_name == CNV_CALL_INTER:
            _default_config = cmd_config.cnv_call_inter_config
        else:
            print(f"\nERROR: Invalid command name: {cmd_name}!\n")
            return
        
        default_config = deepcopy(_default_config)
        convert_config(default_config)

        return default_config

    def _create_config(self, cmd_name):
        """根据命令名称来创建对应命令的用户级配置文件"""

        user_config = dict()

        default_config = self._get_default_config(cmd_name)
        # print(f"===DEBUG===>> default_config: {default_config}")

        base_configurable = default_config['cmd_arguments']['base']
        # smk_configurable = default_config['cmd_arguments']['snakemake'] # 强制所有参数可配置
        tools_configurable = default_config['cmd_arguments']['tools']

        base_first_key = None # 用于后面添加注释
        default_base = default_config['base']
        for arg_k, arg_v in default_base.items():
            if base_configurable == 'all' or arg_k in base_configurable:
                if base_first_key is None:
                    base_first_key = arg_v['arg_name']
                user_config[arg_v['arg_name']] = arg_v['value']

        default_smk = default_config['snakemake']['params']
        user_config['--snakemake'] = dict()
        for arg_k, arg_v in default_smk.items():
            user_config['--snakemake'][arg_v['arg_name']] = arg_v['value']
        
        tools_first_key = None # 用于后面添加注释
        tools_keys = list()
        default_tools = default_config['tools']
        for arg_k, arg_v in default_tools.items():
            if tools_configurable == 'all' or arg_k in tools_configurable:
                if tools_first_key is None:
                    tools_first_key = arg_v['arg_name']
                user_config[arg_v['arg_name']] = arg_v['params']
                tools_keys.append(arg_v['arg_name'])

        # 添加注释、空行
        user_config = config_loader.convert_to_commented_map(user_config)
        user_config.yaml_set_start_comment(comment=f"The configuration file of command <{cmd_name}>!\n")
        user_config.yaml_set_comment_before_after_key(base_first_key, before="\n\nBasic settings")
        user_config.yaml_set_comment_before_after_key('--snakemake', before="\n\nSettings for snakemake")
        user_config.yaml_set_comment_before_after_key(tools_first_key, before="\n\nSettings for tools")
        tools_keys.remove(tools_first_key)
        for k in tools_keys:
            user_config.yaml_set_comment_before_after_key(key=k, before='\n')

        user_config_path = os.path.join(BASE_DIR, 'config', f"user.{cmd_name}_config.yml")
        config_loader.file_dump(user_config, user_config_path)
        # print(f"===DEBUG===>> user_config: {user_config}")

