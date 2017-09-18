'''This script contains functions to concatenate BLAST/usearch output 
with MySQL database information.'''

import MySQLdb

def sigHits(outfile):
	# Builds a list of blastx hits with e < 10^-5
	cdef list hits = []
	with open(outfile, "r") as bx:
		for line in bx:
			line = line.split("\t")
			if float(line[10]) <= 0.00001:
				hits.append(line[0])
	return hits

def subsetSig(infile, hits, outdir):
	# Subsets all entries from hits into a new fasta
	cdef int keep = 0
	cdef str query
	cdef str line
	cdef float count = 0.0
	cdef int total
	total = len(hits)
	query = outdir + "sigHits.fna"
	print("\tSubsetting significant hits...")
	with open(infile, "r") as fasta:
		with open(query, "w") as output:
			for line in fasta:
				if line[0] == ">":
					if " " in line:
						# Select first item if there is a space in the header
						line = line.split()[0]
					if line[1:].strip() in hits:
						output.write(line)
						keep = 1
						count += 1.0
				elif keep == 1:
					output.write(line.upper())
					keep = 0
	return query

#-----------------------------------------------------------------------------

def extractTable(cursor, table):
	# Returns dict of table entries
	t = {}
	cdef str sql
	if table == "Annotations":
		sql = ("SELECT * FROM {};").format(table)
		cursor.execute(sql)
		results = cursor.fetchall()
		for i in results:
			try:
				t[i[0]] = i[1:]
			except IndexError:
				pass
	elif table == "Genes":
		sql = ("SELECT ProteinID, Product FROM {};").format(table)
		cursor.execute(sql)
		results = cursor.fetchall()
		for i in results:
			try:
				t[i[0]] = i[1]
			except IndexError:
				pass
	return t

def concatenateX(cursor, bx, genomes):
	# Extract from Genes table
	cdef list prot = []
	cdef list out
	cdef int st
	cdef int sp
	cdef str string
	p = extractTable(cursor, "Genes")
	for i in bx.keys():
		# Assemble blastx output and append bx without accession/proteinid
		g = genomes[bx[i][0]]
		prod = p[bx[i][1]]
		out = [i, bx[i][0], bx[i][1], prod, g[0], g[1], g[2], g[3], g[4]+"-"+g[5]]
		out.extend(bx[i][2:-2])
		# Add query coverage (match length/query sequence length) and subset sequence
		out.append(int(bx[i][3])/bx[i][-1])
		out.append(bx[i][-2])
		string = ""
		for i in out:
			string += str(i).replace(",", ":") + ","
		prot.append(string[:-1])
	return prot

def concatenateN(cursor, bn):
	# Extracts data from table and returns combined string
	cdef list dna = []
	cdef list out
	cdef str string
	genomes = extractTable(cursor, "Annotations")
	# Assemble blastn output and append bn without accession
	for i in bn.keys():
		g = genomes[bn[i][0]]
		try:
			seq = g[7]
			gc = (seq.count("G") + seq.count("C"))/int(g[3])
		except IndexError:
			gc = "NA"
		out = [i, bn[i][0], g[0], g[1], g[2], g[3], gc]
		out.extend(bn[i][1:-2])
		# Add query coverage (match length/query sequence length) 
		out.append(int(bn[i][2])/bn[i][-1])
		out.append(bn[i][-2])
		string = ""
		for i in out:
			string += str(i).replace(",", ":") + ","
		dna.append(string[:-1])
	return dna, genomes

def printOutput(db, outdir, blastn, blastx):
	# Concatenates output and prints to file
	cdef str header
	cdef str blastcol
	cdef str nhead
	cdef str xhead
	cdef str i
	cursor = db.cursor()
	blastcol = "%ID,MatchLength,#Mismatches,#GapOpenings,QueryStart,QueryEnd,\
HitStart,HitEnd,e,BitScore,QueryCoverage%,Sequence\n"
	header = "Query,Accession,ID,Name,Organism,Description,Taxonomy,SourceLength(bp),"
	nhead = header.replace("ID,Name,", "") + "%GC," + blastcol
	xhead = header + "MoleculeType," + blastcol
	print("\tConcatenating blastn output...")
	dna, genomes = concatenateN(cursor, blastn)
	with open(outdir + "nucleotideResults.csv", "w") as output:
		# Print blastn output
		output.write(nhead)
		for i in dna:
			output.write(i + "\n")
	print("\tConcatenating singificant blastx output...")
	prot = concatenateX(cursor, blastx, genomes)
	with open(outdir + "significantProteinResults.csv", "w") as output:
		output.write(xhead)
		for i in prot:
			output.write(i + "\n")

def getSeqs(infile, blastn, blastx):
	# Extracts gene sequence from input file and appends to dict
	cdef str n
	cdef str line
	cdef int st
	cdef int sp
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
				if n not in blastx.keys():
					n = ""
			elif n:
				# Append length to each dict
				l = len(line)
				st = int(blastx[n][6]) - 1
				sp = int(blastx[n][7]) - 1
				# Determine orientation of match
				if st < sp:
					seq = line[st:sp+1]
				else:
					seq = line[sp:st+1]
				blastx[n].extend([seq, l])
				if n in blastn.keys():
					# Append to blastn
					st = int(blastn[n][5]) - 1
					sp = int(blastn[n][6]) - 1
					blastn[n].extend([line[st:sp+1], l])
	return blastn, blastx

def readBlast(x, n):
	# Reads blastn and blastx input into dictionaries
	cdef str i
	blastn = {}
	blastx = {}
	print("\n\tReading blast output...")
	with open(n, "r") as res:
		for line in res:
			line = line.strip().split()
			if "-" in line[1]:
				blastn[line[0]] = [line[1].split("-")[0]]
			else:
				blastn[line[0]] = [line[1]]
			for i in line[2:]:
				blastn[line[0]].append(i)
	with open(x, "r") as res:
		for line in res:
			line = line.strip().split()
			if float(line[10]) <= 0.00001:
				# Only selet significant blastx hits
				# Append Accession and ID
				blastx[line[0]] = line[1].split("-")
				for i in line[2:]:
					blastx[line[0]].append(i)
	return blastn, blastx
