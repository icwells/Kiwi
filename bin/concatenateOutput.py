'''This script will convert an NCBI flatfile to fasta and upload the data to 
a MySQL database.'''

import argparse
import MySQLdb
from getpass import getpass
from dbIO import extractRow

def extractN(cursor, acc):
	# Extract from Annotations table
	n = extractRow(cursor, "Annotations", acc, "Accession")
	try:
		n = n[0]
		seq = n[8]
		nl = len(seq)
		gc = (seq.count("G") + seq.count("C"))/nl
	except IndexError:
		nl = "NA"
		gc = "NA"
		n = []
		for i in range(9):
			n.append(nl)
	return n, nl, gc

def concatenateX(cursor, name, bx):
	# Extract from Genes table
	n, nl, gc = extractN(cursor, bx[0])
	x = extractRow(cursor, "Genes", bx[1], "ProteinID")
	try:
		x = x[0]
		xl = len(x[9])
	except IndexError:
		xl = "NA"
		x = []
		for i in range(10):
			x.append(xl)
	# Assemble blastx output and append bx without accession/proteinid
	out = [name, bx[0], bx[1], x[4], n[1], n[2], n[3], xl, n[5]]
	out.extend(bx[2:])
	# Add query coverage and subset sequence
	try:
		out.append(int(bx[3])/xl)
	except TypeError:
		out.append("NA")
	st = int(bx[8]) - 1
	sp = int(bx[9]) - 1
	if st < sp:
		out.append(x[9][st:sp+1])
	else:
		out.append(x[9][sp:st+1])
	string = ""
	for i in out:
		string += str(i) + ","
	return string[:-1]

def concatenateN(cursor, name, bn):
	# Extracts data from table and returns combined string
	n, l, gc = extractN(cursor, bn[0])
	# Assemble blastn output and append bn without accession
	out = [name, bn[0], n[1], n[2], n[3], l, gc]
	out.extend(bn[1:-2])
	# Add query coverage
	out.extend([int(bn[2])/bn[-1], bn[-2]])
	string = ""
	for i in out:
		string += str(i) + ","
	return string[:-1]

def printOutput(db, outdir, blastn, blastx):
	# Concatenates output and prints to file
	cursor = db.cursor()
	blastcol = "%ID,MatchLength,#Mismatches,#GapOpenings,QueryStart,QueryEnd,\
HitStart,HitEnd,e,BitScore,QueryCoverage%,Sequence\n"
	header = "Query,Accession,ID,Name,Organism,Description,Taxonomy,SourceLength,"
	nhead = header.replace("ID,Name,", "") + "%GC," + blastcol
	xhead = header + "MoleculeType," + blastcol
	print("\tConcatenating blastn output...")
	with open(outdir + "blastNResults.csv", "w") as output:
		# Print blastn output
		output.write(nhead)
		for i in blastn.keys():
			string = concatenateN(cursor, i, blastn[i])
			if string:
				output.write(string + "\n")
	print("\tConcatenating singificant blastx output...")
	with open(outdir + "significantBlastXResults.csv", "w") as output:
		output.write(xhead)
		for i in blastx.keys():
			string = concatenateX(cursor, i, blastx[i])
			if string:
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
				st = int(blastn[n][5]) - 1
				sp = int(blastn[n][6]) - 1
				blastn[n].extend([line[st:sp+1], len(line)])
	return blastn

def readBlast(indir):
	# Reads blastn and blastx input into dictionaries
	blastn = {}
	blastx = {}
	print("\n\tReading blast output...")
	with open(indir + "blastn.outfmt6", "r") as res:
		for line in res:
			line = line.strip().split()
			if "-" in line[1]:
				blastn[line[0]] = [line[1].split("-")[0]]
			else:
				blastn[line[0]] = [line[1]]
			for i in line[2:]:
				blastn[line[0]].append(i)
	with open(indir + "blastx.outfmt6", "r") as res:
		for line in res:
			line = line.strip().split()
			if line[10] <= 0.00001:
				# Only selet significant blastx hits
				# Append Accession and ID
				blastx[line[0]] = line[1].split("-")
				for i in line[2:]:
					blastx[line[0]].append(i)
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
