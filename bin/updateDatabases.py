'''This script will download current NCBI refSeqs, upload any new entries to 
the MySQL database, and create new BLAST databases.'''

import argparse
import MySQLdb
import os
from subprocess import Popen
from shlex import split
from datetime import datetime
from getpass import getpass
from dbIO import *
from flatFileClass import *
from blastSeqs import makeDB
from version import version

def findNewAcc(old, outfile):
	# Records new accessions in updated RefSeq flat file
	newids = []
	with open(outfile, "r") as flatfile:
		for line in flatfile:
			newline = re.sub(" +", " ", line.strip())
			if "ACCESSION" in newline:
				# Isolate accession number
				acc = newline.split()[1]
				if acc not in old:
					newids.append(acc)
	return newids

def newEntries(db, outfile):
	# Subset new entries to reduce upload time
	refseq = "/tmp/viral.1.genomic.gbff"
	print("\n\tIdentifying new entries...")
	ids = getAccession(db, "Annotations")
	newids = findNewAcc(ids, refseq)
	subFlatFile(refseq, outfile, newids)
	os.remove(refseq)

def updateData(db, outdir, upload, extract):
	# Updates MySQL database and blast databases
	sub = "/tmp/flatFile.subset.gbff"
	if upload == False:
		newEntries(db, sub)
	if extract == False:
		# Add new entries
		columns = getColumns()
		convertFF(sub, db, "Genes", columns, [])
	# Append new data to fasta database files
	print("\tExtracting DNA sequences from database...")
	extractDNA(db, outdir)
	print("\tExtracting protein sequences from database...")
	extractProtein(db, "Genes", outdir)
	# Recompile blast databases
	makeDB(outdir)
	ublastDb(outdir)

def downloadRefSeqs():
	# Downloads and concatenates NCBI viral ref seqs
	print("\n\tDownloading RefSeqs from NCBI...")
	os.chdir("/tmp")
	path = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/"
	v1 = "viral.1.genomic.gbff.gz"
	v2 = "viral.2.genomic.gbff.gz"
	ftp1 = Popen(split("wget " + path + v1))
	ftp1.communicate()
	ftp2 = Popen(split("wget " + path + v2))
	ftp2.communicate()
	# Unzip files and concatenate
	print("\tUnzipping and concatenating files...")
	uz1 = Popen(split("gzip -d " + v1))
	uz1.communicate()
	uz2 = Popen(split("gzip -d " + v2))
	uz2.communicate()
	v1 = v1.replace(".gz", "")
	v2 = v2.replace(".gz", "")
	with open(v1, "a") as output:
		with open(v2, "r") as refseq:
			for line in refseq:
				output.write(line)
	os.remove(v2)

def main():
	starttime = datetime.now()
	wd = os.getcwd()
	parser = argparse.ArgumentParser(description = "This script will download \
current NCBI refSeqs, upload any new entries to the MySQL database,  and \
create new BLAST databases. It can be resumed at each major step (although you \
must call blastSeqs.py directly to continue with making new blast databases).")
	parser.add_argument("-v", action = "store_true", 
help = "Prints version info and exits.")
	parser.add_argument("--id", action = "store_true",
help = "Resumes script from identification phase (skips download).")
	parser.add_argument("--upload", action = "store_true", help = "Resumes \
script from uploading to MySQL database (requires subset flat file).")
	parser.add_argument("--extract", action = "store_true", help = "Resumes \
script from extracting from MySQL database (new entries must have been uploaded).")
	parser.add_argument("-u", help = "MySQL username.")
	parser.add_argument("-d", help = "Path to blast database directory.")
	args = parser.parse_args()
	if args.v:
		version()
	# Get username
	if args.u:
		username = args.u
	else:
		username = input("\tEnter MySQL username: ")
	# Get password
	password = getpass(prompt = "\tEnter MySQL password: ")
	try:
		# Connect to database
		db = MySQLdb.connect("localhost", username, password, "ASUviralDB")
	except:
		print("\nIncorrect password. Access denied.")
		quit()
	if args.d[-1] != "/":
		args.d += "/"
	# Skip steps previous to one specified
	if args.extract == True:
		args.upload = True
	if args.upload == True:
		args.id = True
	if args.id == False:
		downloadRefSeqs()
		os.chdir(wd)
	updateData(db, args.d, args.upload, args.extract)
	db.close()
	print("\tFinished updating datbases.")
	print(("\tTotal runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
