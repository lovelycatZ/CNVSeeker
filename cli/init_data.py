import os

from cli.utils import YamlLoader, update_user_config

current_dir = os.path.abspath(os.path.dirname(__file__))
base_dir = os.path.realpath(os.path.join(current_dir, '../'))
cache_dir = os.path.join(base_dir, '.cnvseeker')
# print("============= DEBUG ===============\n", f"cache dir: {cache_dir}")
if not os.path.exists(cache_dir):
    os.makedirs(cache_dir)

config_loader = YamlLoader()

config_path = os.path.join(base_dir, 'config/config.yml')


#### 统一获取 global_config 配置文件中的值
global_config = config_loader.file_load(config_path, update_dcit={'BASE_DIR':base_dir})

## 基本路径
BASE_DIR = global_config['BASE_DIR']
CNV_DIR = global_config['CNV_DIR']

## 子命令名称
CNV_CALL = global_config['CNV-CALL']
CNV_INTER = global_config['CNV-INTER']
CNV_CALL_INTER = global_config['CNV-CALL-INTER']

CONF_HELPER = global_config['CONF_HELPER']

class CmdConfig(object):
    
    def __init__(self) -> None:

        ## NOTICE: 所有子命令的配置文件对象，在传递给函数时，属于引用传递，任何一个地方的修改会影响原值！

        self.cnv_call_config = self.get_config('cnv_call_config', CNV_CALL)
        self.cnv_inter_config = self.get_config('cnv_inter_config', CNV_INTER)
        self.cnv_call_inter_config = self.get_config('cnv_call_inter_config', CNV_CALL_INTER)

    def get_config(self, config_name, cmd_name):
        # 读取默认配置（CmdConfig 中的 config 对象，最好保留模板字符串）
        config_path = global_config[config_name]
        config = config_loader.file_load(config_path, convert="manual")
        # 更新用户的配置
        config_user = os.path.join(BASE_DIR, 'config', f"user.{cmd_name}_config.yml")
        if os.path.isfile(config_user):
            user_config = config_loader.file_load(config_user)
            update_user_config(config, user_config)
        return config
    
cmd_config = CmdConfig()
