'''This script will convert an NCBI flatfile to fasta and upload the data to 
a MySQL database.'''

import argparse
import MySQLdb
from datetime  import datetime
from getpass import getpass
from blastResults import *
from version import version

def main():
	starttime = datetime.now()
	parser = argparse.ArgumentParser()
	parser.add_argument("-v", action = "store_true", 
help = "Prints version info and exits.")
	parser.add_argument("-u", help = "MySQL username.")
	parser.add_argument("-i", help = "Path to fasta file of query sequences.")
	parser.add_argument("-o", 
help = "Path to input directory (output directory from blastSeqs).")
	parser.add_argument("--blast", action = "store_true", default = False,
help = "Indicates BLAST was run in the previous step (default is ublast).")
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
	if args.o [-1] != "/":
		args.o += "/"
	if args.blast == True:
		x = args.o + "blastx.outfmt6"
		n = args.o + "blastn.outfmt6"
	else:
		x = args.o + "ublastX.outfmt6"
		n = args.o + "ublastN.outfmt6"
	blastn, blastx = readBlast(x, n)
	blastn, blastx = getSeqs(args.i, blastn, blastx)
	printOutput(db, args.o, blastn, blastx)
	db.close()
	print(("\tFinished in {}\n").format(datetime.now() - starttime))

if __name__ == "__main__":
	main()
