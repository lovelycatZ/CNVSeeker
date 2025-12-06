import re



def convert_config(dict_obj, top=None, key_chain_list=[]):
    """自动替换字典中的模板字符串"""

    if top is None:
        top = dict_obj # top 保存原始的 config 文件

    for key, value in dict_obj.items():
        _key_chain_list = key_chain_list.copy()

        if isinstance(value, dict):
            _key_chain_list.append(key)
            convert_config(value, top, _key_chain_list)
            
        if isinstance(value, str):
            tmpls = re.findall(r'\$\{[\.\w-]+?\}', value)
            for tmpl in tmpls:
                keys = tmpl[2:-1].split(".")
                val = top
                for k in keys:
                    val = val[k]
                if not isinstance(val, str):
                    val = str(val)
                value = value.replace(tmpl, val)

            update_dict = top # 用于 top 的更新
            for k in key_chain_list:
                update_dict = update_dict[k]
            update_dict[key] = value # 替换更新

        # TODO int, float, list, ...
