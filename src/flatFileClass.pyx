'''This script defines a class for parsing NCBI flat files.'''

import os
import re
from string import digits as DIGITS
from datetime import datetime
from dbIO import updateDB

class flatFileEntry():
	# Stores relevant data from flat file entry and performs i/o operations
	def __init__(self, line):
		# line should be first line of flatfile entry
		cdef list splt
		line = line.strip()
		line = re.sub(" +", " ", line)
		splt = line.split()
		self.Accession = splt[1]
		self.Length = splt[2]
		self.Type = splt[4]
		self.Shape = splt[5]
		self.Organism = "NA"
		self.Definition = "NA" 
		self.Hierarchy = ""
		self.DNA = "NA"
		self.Date = str(datetime.now())
		# {geneid: {coordinates, protein id, product, strand, exons, translation}}
		self.genes = {}
		self.geneid = ""
		self.cds = {"coor": "NA", "name": "NA", "pid": "NA", 
					"prod": "NA", "strand": "NA", "exons": "NA", "prot": "NA"}
		self.conditions = {"def": 0, "org": 0, "hier": 0, "cds": 0, "trans": 0}

	def geneDict(self):
		# Resets self.cds with blank dictionary
		self.cds = {"geneid": "NA", "coor": "NA", "name": "NA", "pid": "NA", 
						"prod": "NA", "strand": "NA", "exons": "NA", "prot": "NA"}

	def parseLine(self, line):
		# Extract data from remaining lines
		line = line.strip()
		line = re.sub(" +", " ", line)
		# Identify most common lines first
		if line[0] in DIGITS:
			# Remove numbers and spaces
			seq = re.sub(r"\d+", "", line)
			seq = seq.replace(" ", "")
			if line.split()[0] == "1":
				self.DNA = seq.upper()
			else:
				self.DNA += seq.upper()
		elif self.conditions["trans"] == 1:
			# Append remaining lines of translation
			if "gene" in line or "ORIGIN" in line:
				self.genes[self.geneid] = self.cds
				# Reset for next gene
				self.geneDict()
				self.geneid = ""
				self.conditions["trans"] = 0
			else:
				self.cds["prot"] += line.replace('"', '')
		elif line[0] == "/":
			if "/locus_tag=" in line and self.conditions["cds"] == 1:
				self.cds["name"] = line.split("=")[1].replace('"', '')
			elif "/db_xref=" in line and self.conditions["cds"] == 1:
				self.geneid = line[line.find(":")+1:line.rfind("\"")]
			elif "/product=" in line:
				self.cds["prod"] = line.split("=")[1].replace('"', '').replace(',', '')
			elif "/protein_id=" in line and self.conditions["cds"] == 1:
				self.cds["pid"] = line.split("=")[1].replace('"', '')
			elif "/translation=" in line:
				self.cds["prot"] = line.split('"')[1]
				self.conditions["trans"] = 1
		else:
			# Evaluate one and two line entries
			if "DEFINITION" in line:
				line = line.replace(",", "")
				self.Definition = line.replace("DEFINITION ", "")
				self.conditions["def"] = 1
			elif self.conditions["def"] == 1:
				if "ACCESSION" in line:
					self.conditions["def"] = 0
				else:
					# Append second line of definition
					self.Definition += " " + line.replace(",", "")
			elif "ORGANISM" in line:
				self.Organism = line.replace("ORGANISM ", "")
				self.conditions["org"] = 1
			elif self.conditions["org"] == 1:
				if ";" in line:
					self.conditions["org"] = 0
					self.conditions["hier"] = 1
					line = line.replace(";", ":")
					self.Hierarchy += line.replace(",", ":")
				else:
					self.Organism += line 
			elif self.conditions["hier"] == 1:
				if "REFERENCE" in line:
					self.conditions["hier"] = 0
				else:
					line = line.replace(";", ":")
					self.Hierarchy += line.replace(",", ":")
			elif "CDS" in line:
				# Determine coordinates, strand, and number of exons and reformat
				self.conditions["cds"] = 1
				line = line.split()[-1].replace("(", "")
				line = line.replace(")", "")
				if "complement" in line:
					self.cds["strand"] = "-"
					line = line.replace("complement", "")
				if "join" in line:
					self.cds["exons"] = str(line.count(","))
					line = line.replace("join", "")
					line = line.replace(",", "::")
				self.cds["coor"] = line					

	def dbUpload(self, db, table, columns):
		# Formats output for database and uploads
		cdef int passed = 1
		cdef str annotation = ""
		cdef str sequences = ""
		for i in [self.Accession, self.Organism, self.Definition, self.Hierarchy, 
				self.Length, self.Type, self.Shape, self.Date, self.DNA]:
			annotation += '"' + i + '",'
		# Upload genome information to Annotation table
		uploaded = updateDB(db, "Annotations", columns[0], annotation[:-1])
		if uploaded == 0:
			passed = 0
		for j in self.genes:
			coor = self.genes[j]["coor"]
			name = self.genes[j]["name"]
			pid = self.genes[j]["pid"]
			prod = self.genes[j]["prod"]
			strand = self.genes[j]["strand"]
			exons = self.genes[j]["exons"]
			prot = self.genes[j]["prot"]
			# Upload gene specific data to table
			for i in [self.Accession, name, j, pid, prod, coor, strand,
					 exons, self.Date, prot]:
				sequences += '"' + i + '",'
			uploaded = updateDB(db, table, columns[1], sequences[:-1])
			sequences = ""
			if uploaded == 0:
				# Passed will be false if any uploades fail
				passed = 0
		return passed

#-----------------------------------------------------------------------------

def convertFF(infile, db, table, columns, ids):
	# Builds dict of flat file data
	cdef list splt
	cdef int count = 0
	cdef int total = 0
	cdef int uploaded = 1
	# Initialize log file
	cdef str logfile = infile[:infile.rfind("/")+1] + "failedUploads.txt"
	with open(logfile, "w") as log:
		pass
	# Parse flat file
	print("\n\tUploading flat file to database...")
	with open(infile, "r") as flatfile:
		for line in flatfile:
			if line.strip():
				if line.strip() != "//":
					# Extract data from line and append to class
					if "LOCUS" in line:
						genome = flatFileEntry(line)
					elif genome.Accession not in ids:
						genome.parseLine(line)
				elif line.strip() == "//":
					total += 1
					# Upload to database
					if genome.Accession not in ids and "Bacteria: " not in genome.Hierarchy:
						uploaded = genome.dbUpload(db, table, columns)
						if uploaded == 1:
							# Tally successful
							count += 1
						if uploaded == 0:
							# Record accessions of genes which failed 
							with open(logfile, "a") as log:
								log.write(("{}\n").format(genome.Accession))
						uploaded = 1
	print(("\t{} of {} sequences successfuly uploaded.").format(count, total))

def subFlatFile(infile, outfile, ids):
	# Subsets entries to new file
	cdef int save = 0
	cdef str line
	cdef str newline
	cdef str acc
	print("\tSubsetting flat file...")
	with open(infile, "r") as flatfile:
		with open(outfile, "w") as output:
			for line in flatfile:
				newline = re.sub(" +", " ", line.strip())
				if "LOCUS" in newline:
					# Isolate accession number
					acc = newline.split()[1]
					if acc in ids:
						save = 1
				if save == 1:
					# Write all information from failed entry
					output.write(line)
					if newline == "//":
						save = 0
