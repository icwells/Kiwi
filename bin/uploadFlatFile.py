'''This script will convert an NCBI flatfile to fasta and upload the data to 
a MySQL database.'''

import argparse
import MySQLdb
from datetime import datetime
from getpass import getpass
from dbIO import *
from flatFileClass import *

def main():
	starttime = datetime.now()
	parser = argparse.ArgumentParser()
	parser.add_argument("-u", help = "MySQL username.")
	parser.add_argument("--new", action = "store_true",
 help = "Initializes a new table.")
	parser.add_argument("-t", help = "Name of table to upload data to")
	parser.add_argument("-i", help = "Path to input flat file.")
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
	if args.new:
		# Create new table
		newTable(db, args.t)
		ids = []
	else:
		# Get ids from existing table
		ids = getAccession(db, args.t)
	columns = getColumns()
	convertFF(args.i, db, args.t, columns, ids)
	db.close()
	print(("\tTotal runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
