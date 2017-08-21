'''This script will convert an NCBI flatfile to fasta and upload the data to 
a MySQL database.'''

import argparse
import MySQLdb
from getpass import getpass
from dbIO import extractRow

def printOutput(db, outdir, blastn, blastx):
	# Concatenates output and prints to file
	# Make header
	cursor = db.cursor()
	blastcol = "%ID,BitScore,e,HitStart,HitEnd,QueryStart,QueryEnd,Sequence,SeqLength"
	header = "Query,Accession,Gene,Product,Description,Organism,Taxonomy,Length(bp),GC%,"
	header += "GeneID," + blastcol
	header += ",ProteinID," + blastcol + "\n"
	with open(outdir + "significantBlastResults.csv", "w") as output:
		output.write(header)
		for i in blastn.keys():
			# Extract data from each table
			n = extractRow(cursor, "Annotations", blastn[i][0], "Accession")
			n = n[0]
			x = extractRow(cursor, "Genes", blastx[i][0], "ProteinID")
			x = x[0]
			# Calculate length, gc content, and build list of results
			l = len(n[8])
			gc = (n[8].count("G") + n[8].count("C"))/l
			outlist = [i, blastn[i][0], x[1], x[4], n[2], n[1], n[3], l, gc, x[2]]
			# Append all from blastn except query name
			for j in blastn[i][1:]:
				outlist.append(j)
			# Append protein id, all from blastx, protein seq, and length
			outlist.append(x[3])
			for j in blastx[i][1:]:
				outlist.append(j)
			outlist.extend([x[9], len(x[9])])
			# Build string and write
			string = outlist[0]
			for j in outlist[1:]:
				string += "," + str(j)
			output.write(string + "\n")

def getSeqs(infile, blastn):
	# Extracts gene sequence from input file and appends to dict
	n = ""
	with open(infile, "r") as fasta:
		for line in fasta:
			line = line.strip()
			if line[0] == ">":
				# Isolate name and ensure it is in blast dict
				if " " in line:
					n = line.split()[0][1:]
				else:
					n = line[1:].strip()
				if n not in blastn.keys():
					n = ""
			elif n:
				st = int(blastn[n][6])
				sp = int(blastn[n][7])
				seq = line[st:sp+1]
				blastn[n].extend([seq, len(seq)])
	return blastn
			

def readBlast(indir):
	# Reads blastn and blastx input into dictionaries
	blastn = {}
	blastx = {}
	with open(indir + "blastn.outfmt6", "r") as res:
		for line in res:
			line = line.split()
			# query = [Accession, %ID, bit score, e, hit start, hit end, query start, query end]
			blastn[line[0]] = [line[1].split("-")[0], line[2], line[11].strip(),
							line[10], line[8], line[9], line[6], line[7]]
	with open(indir + "blastx.outfmt6", "r") as res:
		for line in res:
			line = line.split()
			if line[0] in blastn.keys():
				# Only selet significant blastx hits
				# query = [protien id, %ID, bit score, e, hit start, hit end, query start, query end]
				blastx[line[0]] = [line[1].split("-")[1], line[2], line[11].strip(), 
								line[10], line[8], line[9], line[6], line[7]]	
	return blastn, blastx

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-u", help = "MySQL username.")
	parser.add_argument("-i", help = "Path to fasta file of query sequences.")
	parser.add_argument("-o", 
help = "Path to input directory (output directory from blastSeqs).")
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
	if args.o [-1] != "/":
		args.o += "/"
	blastn, blastx = readBlast(args.o)
	blastn = getSeqs(args.i, blastn)
	printOutput(db, args.o, blastn, blastx)
	db.close()

if __name__ == "__main__":
	main()
