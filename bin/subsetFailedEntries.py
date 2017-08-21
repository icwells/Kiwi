'''This script subset failed entries from a previous convertFlatFile run to 
more easily identify the formatting error.'''

import argparse
import re

def getAccessions(outdir):
	# Extracts accessions from failedUploads
	ids = []
	with open(outdir + "failedUploads.txt", "r") as log:
		for line in log:
			line = line.strip()
			if line not in ids:
				ids.append(line)
	return ids	

def subFailed(infile, outdir, ids):
	# Subsets failed entries to new file
	save = False
	ext = infile[infile.rfind(".")+1:]
	with open(infile, "r") as flatfile:
		with open(outdir + "failedEntries." + ext, "w") as output:
			for line in flatfile:
				newline = re.sub(" +", " ", line.strip())
				if "LOCUS" in newline:
					# Isolate accession number
					acc = newline.split()[1]
					if acc in ids:
						save = True
				if save == True:
					# Write all information from failed entry
					output.write(line)
					if newline == "//":
						save = False

def main():
	parser = argparse.ArgumentParser(description = "This script subset \
failed entries from a previous convertFlatFile run to more easily identify \
the formatting error.")
	parser.add_argument("-i", help = "Path to input flat file.")
	parser.add_argument("-o", help = "Path to output directory from previous \
convertFlatFile run. Be sure that failedUploads.txt is present.")
	args = parser.parse_args()
	outdir = args.o
	if outdir[-1] != "/":
		outdir += "/"
	ids = getAccessions(outdir)
	subFailed(args.i, outdir, ids)

if __name__ == "__main__":
	main()
