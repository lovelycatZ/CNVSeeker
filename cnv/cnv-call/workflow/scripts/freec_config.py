import argparse

parser = argparse.ArgumentParser(description='Wrapper for making config of ControlFREEC.')
parser.add_argument('--outputDir', action='store')
parser.add_argument('--chrFiles', action='store')
parser.add_argument('--chrLenFile', action='store')
parser.add_argument('--window', action='store')
parser.add_argument('--ploidy', action='store')
parser.add_argument('--sex', action='store')
parser.add_argument('--breakPointThreshold', action='store')
parser.add_argument('--breakPointType', action='store')
parser.add_argument('--maxThreads', action='store')
parser.add_argument('--sambamba', action='store')
parser.add_argument('--SambambaThreads', action='store')
parser.add_argument('--noisyData', action='store')
parser.add_argument('--printNA', action='store')
parser.add_argument('--input_sample', action='store')
parser.add_argument('--input_sample_format', action='store')
parser.add_argument('--input_sample_orientation', action='store')
parser.add_argument('--out_config', action='store')
args = parser.parse_args()


configfile = '{path}'.format(path = args.out_config)
with open(configfile, "w") as c:
	print("[general]\n", file = c)
	for arg, value in vars(args).items():
		if arg.startswith("input") or arg == "out_config":
			continue
		if arg == "captureRegions":
			print("\n[target]\n", file = c)
		print(' = '.join([arg, value]), file = c)
	print("\n[sample]\n", file = c)
	print("mateFile = {}".format(args.input_sample), file = c)
	print("inputFormat = {}".format(args.input_sample_format), file = c)
	print("mateOrientation = {}".format(args.input_sample_orientation), file = c)
