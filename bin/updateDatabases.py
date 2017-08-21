'''This script will download current NCBI refSeqs, upload any new entries to 
the MySQL database, append new sequences to the appropriate fasta references, 
and regenerate blast databases.'''

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

def formatIDs(ids):
	# Returns a list of a dict for use with parallel function
	newids = {}
	newids["idx"] = 0
	newids["ids"] = ids
	return [newids]

def updateData(db, username, password, outdir):
	# Updates MySQL database and blast databases
	ids = getAccession(db, "Annotations")
	columns = getColumns()
	# Add new entries using Accession
	convertFF("/tmp/viral.1.genomic.gbff", db, "Genes", columns, ids)
	# Append new data to fasta database files
	print("\n\tExtracting DNA sequences from database...")
	acc = getAccession(username, password, "Annotations", new = True)
	acc = formatIDs(acc)
	extractDNA(db, outdir, acc)
	print("\n\tExtracting protein sequences from database...")
	pids = getAccession(db, "Genes", "ProteinID", new = True)
	pids = formatIDs(pids)
	extractProtein(username, password, "Genes", outdir, pids)
	# Recompile blast databases
	fna = outdir + "viralRefSeq.fna"
	print("\tConstructing BLAST nucleotide datbase...")
	makeDB(fna, "nucl")
	faa = outdir + "viralRefProt.faa"
	print("\tConstructing BLAST protein datbase...")
	makeDB(faa, "prot")

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
	print("\n\tUnzipping and concatenating files...")
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
	parser = argparse.ArgumentParser()
	parser.add_argument("--resume", action = "store_true",
help = "Resumes upload if a previous attempt is interrupted.")
	parser.add_argument("-u", help = "MySQL username.")
	parser.add_argument("-d", help = "Path to blast database directory.")
	args = parser.parse_args()
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
	if args.resume == False:
		downloadRefSeqs()
		os.chdir(wd)
	updateData(db, args.d)
	db.close()
	print("\tFinished updating datbases.")
	print(("\tTotal runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
