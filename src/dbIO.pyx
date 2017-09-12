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
				if t == 0:
					acolumns += line.split()[0] + ","
				elif t == 1:
					if types == True:
						scolumns += line.strip() + ","
					else:
						scolumns += line.split()[0] + ","
	return [acolumns[:-1], scolumns[:-1]]

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

def getPID(cursor, acc):
	# Returns a list of protein ids corresponding to list of accessions
	cdef str sql
	cdef str i
	cdef list ids = []
	print("\tExtracting ProteinIDs...")
	sql = "SELECT Accession, ProteinID FROM Genes;"
	cursor.execute(sql)
	results = cursor.fetchall() 
	for row in results:
		try:
			if row[0] in acc:
				ids.append(row[1])
		except IndexError:
			pass
	return ids

def getBacAcc(cursor):
	# Stores a list of bacterial accessions to file
	cdef str sql
	cdef list acc = []
	sql = "SELECT Accession FROM Annotations WHERE Hierarchy LIKE 'Bacteria%';"
	cursor.execute(sql)
	results = cursor.fetchall()
	for i in results:
		try:
			acc.append(i[0])
		except IndexError:
			pass
	return acc

def extractDNA(db, outdir):
	# Extracts dna sequences and accessions
	cursor = db.cursor()
	cdef str outfile
	cdef str sql
	acc = getBacAcc(cursor)
	print("\tExtracting nucleotide sequences in fasta format...")
	outfile = outdir + "viralRefSeq.fna"
	sql = 'SELECT Accession, DNA FROM Annotations;'
	# Execute the SQL command
	cursor.execute(sql)
	# Fetch row from table
	results = cursor.fetchall()
	with open(outfile, "w") as fasta:
		for row in results:
			try:
				if len(row[1]) > 2 and row[0] not in acc:
					# Skip entries with missing data
					fasta.write((">{}\n{}\n").format(row[0], row[1]))
					acc.append(row[0])
			except IndexError:
				pass

def extractProtein(db, table, outdir):
	# Extracts fasta protein sequences, protein IDs, and accessions
	cursor = db.cursor()
	cdef str outfile
	cdef str sql
	acc = getBacAcc(cursor)
	ids = getPID(cursor, acc)
	print("\tExtracting protein sequences in fasta format...")
	outfile = outdir + "viralRefProt.faa"
	sql = ('SELECT Accession, ProteinID, Protein FROM {};').format(table)
	# Execute the SQL command
	cursor.execute(sql)
	# Fetch row from table
	results = cursor.fetchall()
	with open(outfile, "w") as fasta:
		for row in results:
			try:
				
				if len(row[2]) > 2 and row[1] not in ids:
					# Skip entries with missing data
					fasta.write((">{}-{}\n{}\n").format(row[0], row[1], row[2]))
					ids.append(row[1])
			except IndexError:
				pass
