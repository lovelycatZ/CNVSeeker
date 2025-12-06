from cli.init_data import CNV_CALL, CNV_INTER, CNV_CALL_INTER
from cli.args.argument import Argument


create = Argument(
    "--create",
    "-c",
    type = str,
    choices = [CNV_CALL, CNV_INTER, CNV_CALL_INTER],
    help = "Create default config file."
)

reset = Argument(
    "--reset",
    "-r",
    type = str,
    choices = [CNV_CALL, CNV_INTER, CNV_CALL_INTER],
    help = "Reset the config file."
)

config_helper_args = [
    create, reset
]
