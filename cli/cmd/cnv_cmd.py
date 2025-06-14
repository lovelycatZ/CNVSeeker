import os

from cli.cmd.command import Command
from cli.init_data import cmd_config, BASE_DIR, CNV_DIR, CNV_CALL, CNV_INTER, CNV_CALL_INTER
from cli.args.cnv_args import cnv_call_args_keychain_map, cnv_inter_args_keychain_map, cnv_call_inter_args_keychain_map


class CNVCallCmd(Command):
    """ CNV-CALL command """

    name = CNV_CALL

    def __init__(self, arguments) -> None:
        super().__init__(arguments)

        self.name = CNV_CALL
        self.work_dir = os.path.join(CNV_DIR, self.name)

    def append_arguments(self, parser, arg):
        super().append_arguments(parser, arg)

    def handle(self, args) -> None:
        self.config_path = os.path.join(BASE_DIR, f'.cnvseeker/.{self.name}_config.yaml')

        super().handle(self.name, args, cnv_call_args_keychain_map, cmd_config.cnv_call_config, self.config_path, self.work_dir)


class CNVInterCmd(Command):
    """ CNV-INTER command """

    name = CNV_INTER

    def __init__(self, arguments) -> None:
        super().__init__(arguments)

        self.name = CNV_INTER
        self.work_dir = os.path.join(CNV_DIR, self.name)

    def append_arguments(self, parser, arg):
        super().append_arguments(parser, arg)

    def handle(self, args) -> None:
        self.config_path = os.path.join(BASE_DIR, f'.cnvseeker/.{self.name}_config.yaml')

        super().handle(self.name, args, cnv_inter_args_keychain_map, cmd_config.cnv_inter_config, self.config_path, self.work_dir)


class CNVCallInterCmd(Command):
    """ CNV-CALL-INTER command """

    name = CNV_CALL_INTER

    def __init__(self, arguments) -> None:
        super().__init__(arguments)

        self.name = CNV_CALL_INTER
        self.work_dir = os.path.join(CNV_DIR, self.name)

    def append_arguments(self, parser, arg):
        super().append_arguments(parser, arg)

    def handle(self, args) -> None:
        self.config_path = os.path.join(BASE_DIR, f'.cnvseeker/.{self.name}_config.yaml')

        super().handle(self.name, args, cnv_call_inter_args_keychain_map, cmd_config.cnv_call_inter_config, self.config_path, self.work_dir)
