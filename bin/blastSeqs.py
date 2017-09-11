'''This script will automate running blastx and blastn for putative 
viral sequences.'''

import argparse
import os
from datetime import datetime
from subprocess import Popen
from shlex import split
from multiprocessing import cpu_count
from blastResults import sigHits, subsetSig
from version import version

def ublast(infile, outdir, dbdir, threads):
	# Calls ublast against protein db for all input and against dna for 
	# significant hits
	cont = True
	print("\n\tRunning ublast against protein sequences...\n")
	cmd = ("./usearch -ublast {} -db {}viralRefProt.udb -evalue 1e-5 \
-maxaccepts 1 -maxrejects 5 -threads {} -blast6out {}ublastX.outfmt6").format(infile,
															 dbdir, threads, outdir)
	try:
		us = Popen(split(cmd))
		us.communicate()
		if us.returncode == 0:
			print("\n\tFinished running ublast against protein sequences.")
	except:
		print(("\tError: ublast protein search failure on {}.").format(infile))
		cont = False
	if cont == True:
		hits = sigHits(outdir + "ublastX.outfmt6")
		query = subsetSig(infile, hits, outdir)
		print("\n\tRunning ublast against nucleotide database...\n")
		cmd = ("./usearch -ublast {} -db {}viralRefSeq.udb -evalue 1e-5 \
-maxaccepts 1  -maxrejects 5 -threads {} -strand both -blast6out \
{}ublastN.outfmt6").format(query, dbdir, threads, outdir)
		try:
			us = Popen(split(cmd))
			us.communicate()
			if us.returncode == 0:
				print("\n\tFinished running ublast against DNA sequences.")
		except:
			print(("\tError: ublast DNA search failure on {}.").format(infile))
	return True

def blastSeqs(infile, outdir, dbdir, threads):
	# Calls blastx on all input and blastn on all hits with e < 10^-5
	starttime = datetime.now()
	print("\n\tRunning blastx against protein database...")
	cont = True
	cmd = ("blastx -query {} -db {}viralRefProt.faa -num_threads {} -max_target_seqs 1 \
-outfmt 6 -evalue 0.00001 -out {}blastx.outfmt6").format(infile, dbdir, threads, outdir)
	try:
		bs = Popen(split(cmd))
		bs.communicate()
		if bs.returncode == 0:
			print(("\tBlastx runtime: {}").format(datetime.now()-starttime))
	except:
		print(("\tError: Blastx failure on {}.").format(infile))
		cont = False
	if cont == True:
		hits = sigHits(outdir + "blastx.outfmt6")
		query = subsetSig(infile, hits, outdir)
		starttime = datetime.now()
		print("\n\tRunning blastn against nucleotide database...")
		cmd = ("blastn -query {} -db {}viralRefSeq.fna -num_threads {} -max_target_seqs 1 \
-outfmt 6 -out {}blastn.outfmt6").format(query, dbdir, threads, outdir)
		try:
			bs = Popen(split(cmd))
			bs.communicate()
			if bs.returncode == 0:
				print(("\tBlastn runtime: {}").format(datetime.now()-starttime))
		except:
			print(("\tError: Blastn failure on {}.").format(infile))
	return True

def ublastDb(indir):
	# Makes ublast databases
	print("\tBuilding ublast databases...\n")
	cmd = "./usearch -makeudb_ublast {} -output {}"
	infile = indir + "viralRefSeq.fna"
	outfile = indir + "viralRefSeq.udb"
	try:
		ndb = Popen(split(cmd.format(infile, outfile)))
		ndb.communicate()
	except:
		print("\tError: Failed to build ublast DNA database.")
	infile = indir + "viralRefProt.faa"
	outfile = indir + "viralRefProt.udb"
	try:
		pdb = Popen(split(cmd.format(infile, outfile)))
		pdb.communicate()
	except:
		print("\tError: Failed to build ublast protein database.")

def makeDB(indir):
	# Makes protein and nucleotide blast databases
	infile = indir + "viralRefSeq.fna"
	print("\tConstructing BLAST databases...")
	cmd = ("makeblastdb -in {} -parse_seqids -dbtype {}")
	try:
		mdb = Popen(split(cmd.format(infile, "nucl")))
		mdb.communicate()
	except:
		print("\tError: Failed to build BLAST nucleotide database.")
	infile = indir + "viralRefProt.faa"
	try:
		mdb = Popen(split(cmd.format(infile, "prot")))
		mdb.communicate()
	except:
		print("\tError: Failed to build BLAST protein database.")

def main():
	starttime = datetime.now()
	parser = argparse.ArgumentParser(description = "This script will \
automate running blastx and blastn for putative viral sequences. \
Be sure to export the path to blast on your machine before using.")
	parser.add_argument("-v", action = "store_true", 
help = "Prints version info and exits.")
	parser.add_argument("-i", help = "Path to fasta file of query sequences.")
	parser.add_argument("-o", help = "Path to working/output directory.")
	parser.add_argument("-d", help = "Path to directory containing database files.")
	parser.add_argument("--blast", action = "store_true", default = False,
help = "Runs BLAST (default is ublast).")
	parser.add_argument("-p", type = int, default = 1,
help = "Number of threads to run BLAST or ublast with.")
	parser.add_argument("--ublastdb", action = "store_true", 
help = "Creates new ublast DNA and protein databases. Place source fasta files \
into the directory specified with -d.")
	parser.add_argument("--blastdb", action = "store_true", 
help = "Creates new blast DNA and protein databases. Place source fasta files \
into the directory specified with -d.")
	args = parser.parse_args()
	if args.v:
		version()
	# Add trailing slash to database directory
	if args.d[-1] != "/":
		args.d += "/"
	if args.ublastdb:
		ublastDb(args.d)
	elif args.blastdb:
		makeDB(args.d)
	else:
		if args.o[-1] != "/":
			# Add trailing slash
			args.o += "/"
		if not os.path.isdir(args.o):
			os.mkdir(args.o)
		if args.p > cpu_count():
			# Prevent too many threads from being called
			args.p = cpu_count()
		if args.blast == True:
			done = blastSeqs(args.i, args.o, args.d, args.p)
		else:
			done = ublast(args.i, args.o, args.d, args.p)
		print(("\tTotal runtime: {}\n").format(datetime.now() - starttime))

if __name__ == "__main__":
	main()
