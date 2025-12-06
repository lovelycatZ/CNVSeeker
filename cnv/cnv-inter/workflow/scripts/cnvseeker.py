######################################################################################################################
import argparse

__version__ = '4.0.0'
print(__file__ )
print('Current Version: ', __version__)

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required = True,
					help = 'Input file in BED or VCF format or a CNV record. If BED, the first four columns should be chromosome, start position, '
						 'end position, CNV type (DEL or DUP) and are split by TAB. If CNV record, '
						 'it should be formated like chr1:100-200-DEL or X:200-400-dup.')
# parser.add_argument('-c', '--cores', default = 1, 
# 					help = 'Maximum number of threads to use. Default: 1')
parser.add_argument('-b', '--build', default = "hg38", 
					help = 'Reference genome version')
parser.add_argument('-o', '--output', required = True,
					help = 'Output file in tsv format.')
args = parser.parse_args()
######################################################################################################################

from interpretation.interpret import interpret


if __name__ == "__main__":
	input = args.input
	# cores = args.cores
	build = args.build
	out_file = args.output
	interpret(input, build, out_file)
