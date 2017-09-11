'''This script subset failed entries from a previous convertFlatFile run to 
more easily identify the formatting error.'''

import argparse
from flatFileClass import subFlatFile
from version import version

def getAccessions(outdir):
	# Extracts accessions from failedUploads
	ids = []
	with open(outdir + "failedUploads.txt", "r") as log:
		for line in log:
			line = line.strip()
			if line not in ids:
				ids.append(line)
	return ids

def main():
	parser = argparse.ArgumentParser(description = "This script subset \
failed entries from a previous convertFlatFile run to more easily identify \
the formatting error.")
	parser.add_argument("-v", action = "store_true", 
help = "Prints version info and exits.")
	parser.add_argument("-i", help = "Path to input flat file.")
	args = parser.parse_args()
	if args.v:
		version()
	# Construct outfile
	outdir = args.i[:args.i.rfind("/")+1]
	ext = infile[infile.rfind(".")+1:]
	outfile = outdir + "flatFile.subset." + ext
	ids = getAccessions(outdir)
	subFlatFile(args.i, outfile, ids)

if __name__ == "__main__":
	main()
