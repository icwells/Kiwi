'''This script will automate running blastx and blastn for putative 
viral sequences.'''

import argparse
import os
from subprocess import Popen
from shlex import split
from multiprocessing import cpu_count

def sigHits(outdir):
	# Builds a list of blastx hits with e < 10^-5
	hits = []
	with open(outdir + "blastx.outfmt6", "r") as bx:
		for line in bx:
			line = line.split("\t")
			if float(line[10]) <= 0.00001:
				hits.append(line[0])
	return hits

def subsetSig(infile, hits, outdir):
	# Subsets all entries from hits into a new fasta
	keep = False
	query = outdir + "sigHits.fna"
	with open(infile, "r") as fasta:
		with open(query, "w") as output:
			for line in fasta:
				if line[0] == ">":
					if line[1:].strip() in hits or line.split()[0].strip()[1:] in hits:
						output.write(line)
						keep = True
				elif keep == True:
					output.write(line)
					keep = False
	return query

def blastSeqs(infile, outdir, dbdir, threads):
	# Calls blastx on all input and blastn on all hits with e < 10^-5
	print("\n\tRunning blastx against protein database...")
	cmd = ("blastx -query {} -db {}viralRefProt.faa -num_threads {} -max_target_seqs 1 \
-outfmt 6 -evalue 0.00001 -out {}blastx.outfmt6").format(infile, dbdir, threads, outdir)
	try:
		bs = Popen(split(cmd))
		bs.communicate()
		if bs.returncode == 0:
			print("\tFinished running blastx.")
	except:
		print("\tError: Blastx failure. Exiting")
		quit()
	hits = sigHits(outdir)
	query = subsetSig(infile, hits, outdir)
	print("\n\tRunning blastn against nucleotide database...")
	cmd = ("blastn -query {} -db {}viralRefSeq.fna -num_threads {} -max_target_seqs 1 \
-outfmt 6 -out {}blastn.outfmt6").format(query, dbdir, threads, outdir)
	try:
		bs = Popen(split(cmd))
		bs.communicate()
		if bs.returncode == 0:
			print("\tFinished running blastn.")
	except:
		print("\tError: Blastn failure. Exiting")
		quit()

def makeDB(infile, typ):
	# Makes protein and nucleotide blast databases
	cmd = ("makeblastdb -in {} -parse_seqids -dbtype {}").format(infile, typ)
	try:
		mdb = Popen(split(cmd))
		mdb.communicate()
		if mdb.returncode == 0:
			print("\n\tFinished constructing BLAST database.")
	except:
		print("\tError: Failed to build BLAST database.")

def main():
	parser = argparse.ArgumentParser(description = "This script will \
automate running blastx and blastn for putative viral sequences. \
Be sure to export the path to blast on your machine before using.")
	parser.add_argument("-i", help = "Path to fasta file of query sequences.")
	parser.add_argument("-o", help = "Path to working/output directory.")
	parser.add_argument("-d", 
help = "Path to directory containing blast databases.")
	parser.add_argument("-p", type = int, default = 1,
help = "Number of threads to run BLAST with.")
	parser.add_argument("--nucdb", action = "store_true", help = "Creates a new \
dna blast database. Place source fasta file into the directory specified with -d.")
	parser.add_argument("--protdb", action = "store_true",
help = "Creates a new protein blast database. \
Place source fasta file into the directory specified with -d.")
	args = parser.parse_args()
	# Add trailing slash to database directory
	if args.d[-1] != "/":
		args.d += "/"
	if args.nucdb or args.protdb:
		if args.nucdb:
			infile = args.d + "viralRefSeq.fna"
			print("\tConstructing BLAST nucleotide datbase...")
			makeDB(infile, "nucl")
		if args.protdb:
			infile = args.d + "viralRefProt.faa"
			print("\tConstructing BLAST protein datbase...")
			makeDB(infile, "prot")
	else:
		if args.o[-1] != "/":
			# Add trailing slash
			args.o += "/"
		if args.p > cpu_count():
			# Prevent too many threads from being called
			args.p = cpu_count()
		blastSeqs(args.i, args.o, args.d, args.p)

if __name__ == "__main__":
	main()
