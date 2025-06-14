import re

# from string import Template
from ruamel import yaml


class YamlLoader:
    def __init__(self):
        self.yaml = yaml.YAML()
        # yaml 对象设置
        self.yaml.preserve_quotes = True
        # 设置缩减，mapping 是普通的值的缩进；sequence 是列表的缩进（列表值的首字符，不包括 - ）；offset 是 "-" 的缩进！
        self.yaml.indent(mapping=2, sequence=4, offset=2)

    def file_load(self, file, convert: str = "auto", update_dcit: dict = None):
        """加载配置文件

        prams:
            file: 配置文件路径
            convert: str, 设置是否在读取文件时自动进行模板字符串的转换。可选 {"auto", "manual"}。default: "auto"
            update_dcit: dict, 设置需要进行更新的配置。default: None

        return: 配置信息对象（字典类型的子类）
        """
        with open(file, 'r', encoding='utf-8') as f:
            data = f.read()

        yaml_info = self.yaml.load(data)

        if update_dcit:
            for k,v in update_dcit.items():
                if k in yaml_info:
                    yaml_info[k] = v

        if convert == "auto":
            convert_config(yaml_info)

        return yaml_info

    def file_dump(self, data, new_file):
        with open(new_file, 'w', encoding='utf-8') as f:
            self.yaml.dump(data, f)

    def convert_to_commented_map(self, yaml_info):
        # CommentedMap 是 ruamel 默认的对象
        return yaml.CommentedMap(yaml_info)

def convert_config(dict_obj, top=None, key_chain_list=[]):
    """自动替换字典中的模板字符串"""

    if top is None:
        top = dict_obj # top 保存原始的 config 文件

    for key, value in dict_obj.items():
        _key_chain_list = key_chain_list.copy() # 使用复制的列表进行后续处理，保证原列表由于引用传递不被修改

        if isinstance(value, dict):
            _key_chain_list.append(key)
            convert_config(value, top, _key_chain_list)
            
        if isinstance(value, str):
            # print('======= keys: ', key_chain_list, key)      # === DEBUG ===
            # print('old value: ', value)                       # === DEBUG ===
            tmpls = re.findall(r'\$\{[\.\w-]+?\}', value)
            for tmpl in tmpls:
                keys = tmpl[2:-1].split(".")
                val = top
                for k in keys:
                    val = val[k]
                if not isinstance(val, str):
                    # print(f"(convert_config) Warning: The type of '{tmpl.strip('${} ')}' is not str, but will be forced to the str!")
                    val = str(val)
                value = value.replace(tmpl, val)
            # print('new value: ', value)                       # === DEBUG ===

            update_dict = top # 用于 top 的更新
            for k in key_chain_list:
                update_dict = update_dict[k]
            # print('------- before：', '\n', top)              # === DEBUG ===
            update_dict[key] = value # 替换更新
            # print('------- after：', '\n', top)               # === DEBUG ===

        # TODO int, float, list, ...


def arg_key_convert(arg_key:str, mode:str='var', short_name=False):
    """转换命令行参数的参数名
    prams:
        arg_key: str, 参数的 key
        mode: str, 转换模式，可选 {'var', 'cmd', 'rm_prefix', 'add_prefix'}。default: var
            'var': 把命令行参数中的值转换为合法的变量值：去除前面的 - , 将中间的 - 换成 _;
            'cmd': 和 'var' 相反，把变量值转为命令行中的参数形式。
            'rm_prefix': 只去除前面的 - 或 -- 
            'add_prefix': 如果输入的 arg_key 不以 -/-- 开头，则在其前面添加 -/--。short_name 参数控制是 - 还是 --。
        short_name: 当 mode 为 'cmd' 时，需要说明当前参数名是否为 short name。default: False

    return: str, 返回转换后的参数的参数名
    """
    if mode == 'var':
        arg_key = arg_key.strip().lstrip('-')
        arg_key = arg_key.replace('-', '_')
    elif mode == 'cmd':
        arg_key = arg_key.strip().replace('_', '-')
        arg_key = arg_key.lstrip('-')
        arg_key = '-'+arg_key if short_name else '--'+arg_key
    elif mode == 'rm_prefix':
        arg_key = arg_key.strip().lstrip('-')
    elif mode == 'add_prefix':
        if not arg_key.startswith('-'):
            arg_key = f'-{arg_key}' if short_name else f'--{arg_key}'
    return arg_key

# def parse_cmd_arguments(args_str:str, is_cvt_keys:bool=True) -> dict:
#     """从命令行字符串中提取命令行参数和对应的值"""

#     pattern = re.compile(r'(--?[\w-]+)\s+([\w/][\w/-]*)?') # TODO 1、参数值没有匹配引号
#     res = pattern.findall(args_str)
#     print("================ DEBUG =========== PARSE CMD : ", res)
#     #### 存在的问题：(20240601：已使用更为简单的策略(下方的 parse_cmd_arguments 方法)，问题 1，2 已解决，问题 3 有待验证。)
#     # 1、匹配出错
#     #    python cnvseeker.py test --home-dir maoxinxin/var-pipe/test --fastp 'min-len 222'（这个写法不合法！！）
#     #    匹配到的结果：[('-len', '222')]
#     # 2、action 的值（bool 值）没有进行处理
#     # 3、参数值没有匹配引号（内部的引号）

#     args_dict = {}
#     if res:
#         for arg in res:
#             k,v = arg
#             k = arg_key_convert(k) if is_cvt_keys else k
#             args_dict[k] = v
#     else:
#         print(f'ERROR: Not matched arguments in "{args_str}"')
#         exit(1)
#     return args_dict

def parse_cmd_arguments(args_str:str, is_cvt_keys:bool=True) -> dict:
    """从命令行字符串中提取命令行参数和对应的值

    命令行字符串格式为: name=a,age=20,male (参数之间用逗号隔开，参数值使用等号赋值);
    选项前面可以加前缀 "-" 或 "--"，若不加，则自动根据选项的长度来自动添加前缀
    不处理参数名内部的 "-" 和 "_" 的问题;
    没有区分参数值中的数字和字符串类型，统一当作字符串处理.
    """
    args_dict = {}
    for arg in args_str.strip().split(","):
        if arg.find("=") == -1:
            k, v = arg, True
        elif len([substr.start() for substr in re.finditer('=', arg)]) > 1:
            k, *v = arg.split("=")
            v = '='.join(v)
            v = v.strip('\'"') # 去除引号
        else:
            k, v = arg.split("=")
            v = v.strip('\'"') # 去除引号
        # k = arg_key_convert(k, mode="rm_prefix") if is_cvt_keys else k # 去除选项前面的 -/--
        # 保留选项前面的 -/-- ，即长、短参数由用户自己决定。若用户没有提供前缀，则根据 key 的长度决定 short/long name
        k = arg_key_convert(k, mode='add_prefix', short_name=len(k)==1)
        # 兼容使用多个重复选项名的情况：如果有多个相同的 k，则转换为列表
        if k in args_dict.keys():
            if isinstance(args_dict[k], list):
                args_dict[k].append(v)
            else:
                args_dict[k] = [args_dict[k]]
                args_dict[k].append(v)
        else:
            args_dict[k] = v

    return args_dict


def reverse_parse_cmd_arguments(param_dict):
    """
    从参数字典中拼接命令行字符串。(parse_cmd_arguments 函数的逆过程)
    """
    param_str = list()
    for k,v in param_dict.items():
        if v is True:
            param_str.append(f'{k}')
        elif v is False:
            continue
        elif isinstance(v, list):
            value_list = [f'{k} "{p}"' for p in v]
            param_str += value_list
        else:
            param_str.append(f'{k} "{v}"')

    param_str = ' '.join(param_str)
    return param_str

def update_user_config(default_config: dict, user_config: dict):
    """使用用户的配置文件更新默认配置"""

    default_base = default_config['base']
    default_smk = default_config['snakemake']['params']
    default_tools = default_config['tools']

    for k,v in user_config.items():
        if isinstance(k,str) and k.startswith('-'):
            if not isinstance(v, dict):
                # base
                for arg_k, arg_v in default_base.items():
                    if k == arg_v.get('arg_name') or k == arg_v.get('arg_short_name'):
                        default_base[arg_k]['value'] = v
                        if 'required' in default_base[arg_k] and default_base[arg_k]['required']:
                            default_base[arg_k]['required'] = False
                        break
            elif k == '--snakemake' or k == '-smk':
                # snakemake
                for user_arg_k, user_arg_v in v.items():
                    for arg_k, arg_v in default_smk.items():
                        if user_arg_k == arg_v.get('arg_name') or user_arg_k == arg_v.get('arg_short_name'):
                            default_smk[arg_k]['value'] = user_arg_v
                            break
            else:
                # tools
                for tool,tool_val in default_tools.items():
                    if k == tool_val.get('arg_name') or k == tool_val.get('arg_short_name'):
                        for user_arg_k, user_arg_v in v.items():
                            tool_val['params'][user_arg_k] = user_arg_v
                        break

        elif k in default_config and k not in ['cmd_arguments', 'base', 'snakemake', 'tools']: # TODO 这个地方，写死不是很好，待后续优化
                
                default_config[k] = v
    
    # print(f"===DEBUG===>> default_config: {default_config}")
    return default_config


if __name__ == "__main__":

    ...