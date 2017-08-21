'''These functions will upload and extract data from the ASUviral database.'''

from sys import stdout
import MySQLdb
import os

def updateDB(db, table, columns, values):
	# Updates existing table
	cdef str sql
	cursor = db.cursor()
	# Conctruct SQL statement
	sql = ("INSERT INTO " + table + "(" + columns 
			+ ") VALUES(" + values + ");")
	try:
		# Insert new row
		cursor.execute(sql)
		db.commit()
		return 1
	except:
		db.rollback()
		err = values.split(",")[0].replace('"', '')
		if "ProteinID" in columns:
			err += ": " + values.split(",")[3].replace('"', '')
		print(("\tThere was an error uploading {}").format(err))
		return 0

def buildValues(gene, columns):
	# Organizes input data
	cdef str values = ""
	cdef list cols = columns.split(",")
	cdef str i
	for i in cols:
		try:
			values += '"' + gene[i] + '",'
		except KeyError:
			# Append dash for missing info
			values += '"-",'
	return values[:-1]

def getAccession(db, table, token = "Accession", new = False):
	# Retrieves list of accession numbers from database
	# Add to set to automatically skip repeats, return as list
	ids = set()
	cursor = db.cursor()
	if new == True:
		sql = ("SELECT {} FROM {} WHERE Date >= CURDATE();"
				).format(token, table)
	else:
		sql = ("SELECT {} FROM {};").format(token, table)
	cursor.execute(sql)
	result = cursor.fetchall()
	for i in result:	
		ids.add(i[0])
	return list(ids)

def getColumns(types = False):
	# Build list of two column statements
	cdef str acolumns = ""
	cdef str scolumns = ""
	cdef str line
	cdef int t
	with open("tableColumns.txt") as col:
		for line in col:
			if line[0] == "#":
				if "Annotations" in line:
					t = 0
				elif "Sequences" in line:
					t = 1
			elif line.strip():
				if t == 0 and types == False:
					acolumns += line.split()[0] + ","
				elif t == 1:
					if types == True:
						scolumns += line.strip() + ","
					else:
						scolumns += line.split()[0] + ","
	if types == False:
		acolumns = acolumns[:-1]
	return [acolumns, scolumns[:-1]]

def newTable(db, table):
	# Initializes new table
	cdef list columns
	cursor = db.cursor()
	columns = getColumns(types = True)
	cmd = ("CREATE TABLE {}({});").format(table, columns[1])
	try:
		cursor.execute(cmd)
		db.commit()
	except:
		db.rollback()
		print("\tError creating new table. Exiting.\n")
		quit()

#-----------------------------------------------------------------------------

def extractRow(cursor, table, name, column):
	# Extract data from table
	cdef str sql
	sql = ('SELECT * FROM {} WHERE {} = "{}";').format(table, column, name)
	# Execute the SQL command
	cursor.execute(sql)
	# Fetch row from table
	results = cursor.fetchall()
	return results

def extractCSV(db, table, outdir, ids, columns):
	# Returns all data associated with an accession
	cursor = db.cursor()
	cdef str token
	cdef str i
	cdef str line
	print(("\n\tExtracting from {}...").format(table))
	token = "ProteinID"
	if table == "Annotations":
		token = "Accession"
	with open(outdir + table + ".csv", "w") as output:
		# Write header
		output.write(columns[:columns.rfind(",")] + "\n")
		# Extract data from table
		for i in ids:
			row = extractRow(cursor, table, i, token)
			# Remove sequences from list
			row = row[0][:-1]
			line = ""
			for j in row:
				line += str(j) + ","
			output.write(line[:-1] + "\n")

def appendFasta(infile):
	# Adds file from parallel extraction into result fasta
	cdef str line
	cdef str outfile
	outfile = infile[:infile.find(".")] + infile[infile.rfind("."):]
	with open(outfile, "a") as output:
		with open(infile, "r") as fasta:
			for line in fasta:
				output.write(line)
	os.remove(infile)

def extractDNA(username, password, outdir, para):
	# Extracts fasta dna sequences, gene IDs, and accessions
	db = MySQLdb.connect("localhost", username, password, "ASUviralDB")
	cursor = db.cursor()
	cdef str i
	cdef list ids 
	cdef str outfile
	outfile = outdir + ("viralRefSeq.{}.fna").format(para["idx"])
	ids = para["ids"]
	with open(outfile, "w") as fasta:
		for i in ids:
			row = extractRow(cursor, "Annotations", i, "Accession")
			try:
				row = row[0]
				if len(row[8]) > 2:
					# Skip entries with missing data
					fasta.write((">{}-{}bp\n{}\n").format(row[0], row[4], row[8]))
			except IndexError:
				pass
	appendFasta(outfile)

def extractProtein(username, password, table, outdir, para):
	# Extracts fasta dna sequences, gene IDs, and accessions
	db = MySQLdb.connect("localhost", username, password, "ASUviralDB")
	cursor = db.cursor()
	cdef str i
	cdef list ids
	cdef str outfile
	outfile = outdir + ("viralRefProt.{}.faa").format(para["idx"])
	ids = para["ids"]
	with open(outfile, "w") as fasta:
		for i in ids:
			row = extractRow(cursor, table, i, "ProteinID")
			try:
				row = row[0]
				if len(row[9]) > 2:
					# Skip entries with missing data
					fasta.write((">{}-{}\n{}\n").format(row[0], row[3], row[9]))
			except IndexError:
				pass
	appendFasta(outfile)
