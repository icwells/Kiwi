'''This script will extract data from the ASUviralDB in various formats.'''

import argparse
import MySQLdb
from datetime import datetime, date
from math import ceil
from multiprocessing import Pool, cpu_count
from functools import partial
import os
from subprocess import Popen
from shlex import split
from getpass import getpass
from dbIO import *

def backup():
	# Backup database to local Linux machine
	timestamp = str(date.today())
	f = open(os.devnull, 'w')
	sys.stdout = f
	try:
		p = Popen(split("mysqldump -u root -p --result-file=sampleIndex." 
					+ timestamp + ".sql 'sampleIndex'"), stdout = f)
		# Wait for completion
		p.communicate()
		# Check for errors
		if(p.returncode != 0):
			raise
		print("Backup done for sampleIndex.")
	except:
		print("Backup failed for sampleIndex.")

def main():
	starttime = datetime.now()
	parser = argparse.ArgumentParser(description = 
"This script will extract data from the ASUviralDB in various formats.")
	parser.add_argument("--backup", action = "store_true",
help = "Backup database to local Linux machine.")
	parser.add_argument("-u", help = "MySQL username.")
	parser.add_argument("-t", help = "Name of table to extract data from \
(not needed for extracting dna sequences).")
	parser.add_argument("--csv", action = "store_true",
help = "Extract all data from a given table to a csv (sequences will be omitted).")
	parser.add_argument("--dna", action = "store_true",
help = "Extract dna sequences from database in fasta format.")
	parser.add_argument("--protein", action = "store_true",
help = "Extract protein sequences from database in fasta format.")
	parser.add_argument("-d", help = "Path to output directory.")
	args = parser.parse_args()
	if args.backup == True:
		# Backup database
		backup()
	else:
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
			print("\n\tIncorrect password. Access denied.")
			quit()
		if args.d:
			# Add trailing / to directory if it necessary
			if args.d[-1] != "/":
				args.d += "/"
			if not os.path.isdir(args.d):
				os.mkdir(args.d)
		if args.csv:
			# Extract csv annotation
			columns = getColumns()
			if args.t == "Annotations":
				ids = getAccession(db, args.t)
				columns = columns[0]
			else:
				ids = getAccession(db, args.t, "ProteinID")
				columns = columns[1]
			extractCSV(db, args.t, args.d, ids, columns)
			print("\tFinished writing data to file.\n")			
		else:
			# Extract protein and/or dna sequences
			if args.dna == True:
				extractDNA(db, args.d)
				print("\tFinished writing DNA sequences to file.\n")			
			if args.protein == True:
				extractProtein(db, args.t, args.d)
				print("\tFinished writing protein sequences to file.\n")			
	print(("\tTotal runtime: {}.").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
