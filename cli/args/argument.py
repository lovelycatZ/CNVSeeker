from copy import deepcopy

from cli.utils import arg_key_convert, convert_config


class Argument:
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def add_to_parser(self, parser):
        parser.add_argument(*self.args, **self.kwargs)

    def add_to_group(self, group):
        group.add_argument(*self.args, **self.kwargs)


# TODO 这两个列表可以提取到配置文件中
# 针对 Argument(argparse 的参数) 类的两种输入
pos_args_keys = ['arg_short_name', 'arg_name']
kw_args_keys = ['action', 'type', 'choices', 'required', 'help']

def get_args_kwargs(args_dict:dict, args_keys:list, kwargs_keys:list):
    """根据 args keys 列表和 kwargs keys 列表，从 args dict 中获取对应的参数"""

    args = []
    kwargs = {}
    for k in args_keys:
        if args_dict.get(k):
            args.append(args_dict.get(k))
    for k in kwargs_keys:
        if args_dict.get(k):
            kwargs[k] = eval(args_dict.get(k)) if k=='type' else args_dict.get(k) # === TODO 可以优化，小心注入问题 ===
    return args, kwargs

# 实现从配置文件中自动配置默认参数
def get_argumets_from_config(_config:dict):

    config = deepcopy(_config)
    # 转换模板字符串
    convert_config(config)

    args = {} 
    args_keychain_map = {}
    args_cats = config['cmd_arguments'].keys()

    for args_cat in args_cats:
        if args_cat == 'snakemake':
            args[args_cat] = []
            p, pval = args_cat, config[args_cat]
            if pval is None:  # 这个是兼容配置文件中 snakemake 那个部分没有参数的情况（同下面 2 个参数类型的 'if args_dict is None' 语句）
                continue
            arg_name, arg_property = get_args_kwargs(pval, pos_args_keys, kw_args_keys)

            # 配置 snakemake 默认值
            smk_params_dict = pval['params']
            if not isinstance(smk_params_dict, dict): # pval['params'] 不是字典类型，说明配置文件有误（或 params 为空），直接 continue
                continue
            smk_value = list()
            for smk_p, smk_pval in smk_params_dict.items():
                if smk_pval['value'] is True:
                    smk_value.append(f"{smk_pval['arg_name'].strip()}")
                elif smk_pval['value'] is False:
                    continue
                else: # TODO 其它 list 等特殊类型待后续优化
                    smk_value.append(f"{smk_pval['arg_name'].strip()}={str(smk_pval['value']).strip()}")

            arg_property['default'] = ','.join(smk_value)

            argument = Argument(*arg_name, **arg_property)
            args[args_cat].append(argument)

        elif args_cat == 'base':
            args[args_cat] = []
            args_dict = config[args_cat]
            args_keychain_map['base'] = {}
            configurable_args = config['cmd_arguments']['base']
            if args_dict is None:
                continue
            for p, pval in args_dict.items():
                # 如果不可配置，则跳过
                if configurable_args!="all" and p not in configurable_args:
                    continue
                args_keychain_map['base'][arg_key_convert(pval.get('arg_name'))] = '.'.join(['base', p])

                arg_name, arg_property = get_args_kwargs(pval, pos_args_keys, kw_args_keys)
                argument = Argument(*arg_name, **arg_property)
                args[args_cat].append(argument)

        elif args_cat == 'tools':
            args[args_cat] = []
            args_dict = config[args_cat]
            args_keychain_map['tools'] = {}
            configurable_args = config['cmd_arguments']['tools']
            if args_dict is None:
                continue
            for p, pval in args_dict.items():
                if configurable_args!="all" and p not in configurable_args:
                    continue
                args_keychain_map['tools'][arg_key_convert(pval.get('arg_name'))] = '.'.join(['tools', p])

                arg_name, arg_property = get_args_kwargs(pval, pos_args_keys, kw_args_keys)
                argument = Argument(*arg_name, **arg_property)
                args[args_cat].append(argument)

        else: ...

    return args, args_keychain_map
