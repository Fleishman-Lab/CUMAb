import argparse
from scripts.modules_graft import graft_sequences, find_interface_locs
from scripts.modules_args import read_config

def main():
	# parse arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("sequence", type=str, help="Sequence of the antibody to humanize")
	parser.add_argument("scheme", type=str, help="Numbering scheme of the antibody")
	parser.add_argument("--screens", default=None, help="Screens to use for humanization")
      
	args = parser.parse_args()
	sequence = args.sequence
	scheme = args.scheme
	screens = args.screens

	# create csv of grafted sequences
	graft_sequences(sequence, scheme, screens)

if __name__ == '__main__':
    main()