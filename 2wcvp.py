#!/usr/bin/env python3

### SAFE AND REQUIRED IMPORTS
import ftfy
import getopt
import lzma
import os
import re
import regex ### seems to have a memory leak, do not use
import sys
from thefuzz import fuzz ### pip3 install thefuzz[speedup]
from thefuzz import process



### PRINT TO STANDARD ERROR
def eprint(*args, **kwargs):
	print(*args, file = sys.stderr, **kwargs)



### USER SETTINGS
settings = {}
settings['approximate'] = 'thefuzz'
settings['column'] = 0
settings['debug'] = False
settings['input'] = ''
settings['lower'] = False
settings['original'] = False
settings['verbose'] = False
settings['wcvp'] = ''



### OTHER SETTINGS
settings['ERRORS'] = (0.025, 0.050, 0.075, 0.100)
settings['minScore'] = 90



### READ OPTIONS
approximate = ('none', 'thefuzz', 'regex')
#
# add tre or fuzzy set (https://pypi.org/project/fuzzyset/)?
#
inputError = 'input file to be analyzed for synonymy (header line assumed; required): -i file.tsv | --input=file.tsv'
wcvpError = 'WCVP synonymy file (tested with v5.0; required): -w wcvp.xz | --wcvp=wcvp.xz'
try:
	arguments, values = getopt.getopt(sys.argv[1:], 'a:c:dhi:lovw:', ['approximate=', 'column=', 'debug', 'help', 'input=', 'lower', 'original', 'verbose', 'wcvp='])
except getopt.error as err:
	eprint(str(err))
	sys.exit(2)
for argument, value in arguments:
	if argument in ('-a', '--approximate') and value in approximate:
		settings['approximate'] = value
	elif argument in ('-c', '--column') and int(value) >= 0:
		settings['column'] = int(value)
	elif argument in ('-d', '--debug'):
		settings['debug'] = True
	elif argument in ('-h', '--help'):
		eprint('\n\nA Python3 script for correcting names to WCVP (https://doi.org/10.1038/s41597-021-00997-6)')
		eprint(f"use approximate token matching if exact token matching fails (optional; default = {settings['approximate']}): -a {'|'.join(approximate)} | --approximate={'|'.join(approximate)}")
		eprint(f"input file column containing scientificName (zero indexed; optional; default = {settings['column']}): -c int | --column=int")
		eprint('print debug info to STDERR (optional): -d | --debug')
		eprint(inputError)
		eprint(f"assume input file is all lowercase (optional; default = {settings['lower']}): -l | --lower")
		eprint(f"output both original and corrected name (optional; default = {settings['original']}): -o | --original")
		eprint(f"verbosely insert a column after scientificName detailing correction type (optional; defauilt = {settings['verbose']}): -v | --verbose")
		eprint(wcvpError + '\n')
		sys.exit(0)
	elif argument in ('-i', '--input'):
		if os.path.isfile(value):
			settings['input'] = value
		else:
			eprint(f"input file does not exist {value}")
			sys.exit(2)
	elif argument in ('-l', '--lower'):
		settings['lower'] = True
	elif argument in ('-o', '--original'):
			settings['original'] = True
	elif argument in ('-v', '--verbose'):
		settings['verbose'] = True
	elif argument in ('-w', '--wcvp'):
		if os.path.isfile(value):
			settings['wcvp'] = value
		else:
			eprint(f"wcvp file does not exist {value}")
			sys.exit(2)



### START/END
if not settings['input']:
	eprint(inputError)
	sys.exit(2)
elif not settings['wcvp']:
	eprint(wcvpError)
	sys.exit(2)
else:
	eprint('started...')
	for key, value in settings.items():
		eprint(f"{key} = {value}")



### REGXP
INTRA = re.compile(' (convar\.|f\.|grex|lusus|microgene|modif\.|monstr\.|mut\.|nothosubsp\.|nothovar\.|proles|provar\.|stirps|subf\.|sublusus|subproles|subso|subsp\.|subvar\.|var\.)$')
NOTHO = re.compile('×')



### FUNCTIONS
def agrep(names, pattern, amount):
	best = ['NO MATCH!']
	count = 0
	matches = {} ### match => score
	maximum = 0
	if settings['approximate'] == 'thefuzz':
		matches = dict(process.extract(pattern, names, scorer = fuzz.ratio))
	elif settings['approximate'] == 'regex':
		previousError = 0
		for error in settings['ERRORS']:
			currentError = int(error*len(pattern))
			if currentError == previousError:
				continue
			previousError = currentError
			expression = regex.compile(f"({re.escape(pattern)}){{e<={currentError}}}")
			for name in names:
				match = regex.search(expression, name)
				if match:
					matches[name] = int(100*(1-(sum(match.fuzzy_counts)/len(pattern))))
			if len(matches):
				break
	if len(matches):
		maximum = max(matches.values())
		if maximum >= settings['minScore']:
			best = [key for key, value in matches.items() if value == maximum]
			count = len(best)
	if amount == 'one':
		return count, maximum, best[0]
	elif amount == 'all':
		return count, maximum, best
	else:
		return 0, maximum, f" use all|one! NOT '{amount}'"

def extractGenus(x):
	y = x.split(' ')
	if re.search(NOTHO, y[0]):
		return f"{y[0]} {y[1]}"
	else:
		return y[0]

def insertColumns(algorithm, columns, new):
	i = 1
	if settings['original'] == True:
		columns.insert(settings['column']+1, new)
		i = 2
	else:
		columns[settings['column']] = new
	if settings['verbose'] == True:
		columns.insert(settings['column']+i, algorithm)

def insertSpecies(ambiguous, inputAuthor, inputName, outputAuthor, outputName, synonymy):
	if settings['lower'] == True:
		inputAuthor = inputAuthor.lower()
		inputName = inputName.lower()
	inputNameAuthor = nameAuthor(inputAuthor, inputName)
	outputNameAuthor = nameAuthor(outputAuthor, outputName)
	if len(outputNameAuthor) == 0: ### for 'synonym' without an accepted name (wtf?)
		outputNameAuthor = inputNameAuthor
	if inputName in synonymy and synonymy[inputName] != outputNameAuthor:
		ambiguous[inputName] = True
	if inputNameAuthor in synonymy and synonymy[inputNameAuthor] != outputNameAuthor:
		ambiguous[inputNameAuthor] = True
	synonymy[inputName] = outputNameAuthor
	synonymy[inputNameAuthor] = outputNameAuthor

def insertInfraspecificX(ambiguous, inputAuthor, inputName, outputAuthor, outputName, synonymy):
	insertSpecies(ambiguous, inputAuthor, inputName, outputAuthor, outputName, synonymy)
	insertSpecies(ambiguous, inputAuthor, ' '.join(inputName.split(' ')[0:3]), outputAuthor, outputName, synonymy)

def insertIntraspecific(ambiguous, inputAuthor, inputName, outputAuthor, outputName, synonymy):
	insertSpecies(ambiguous, inputAuthor, inputName, outputAuthor, outputName, synonymy)
	insertSpecies(ambiguous, inputAuthor, ' '.join(inputName.split(' ')[0:2]), outputAuthor, outputName, synonymy)

def nameAuthor(a, n):
	if len(a) > 0:
		return f"{n} {a}"
	else:
		return n

def tokenize(x):
	y = x.split(' ')
	z = []
	start = 1
	if len(y) > 2:
		for k in range(0, 3):
			if re.search(NOTHO, y[k]):
				start += 1
	for k in reversed(range(start, len(y))):
		t = ' '.join(y[0:k+1])
		if not re.search(INTRA, t): 
			z.append(t)
	return z



### WCVP COLUMNS
KEW_ID = 0
FAMILY = 1
GENUS = 2
SPECIES = 3
INFRASPECIES = 4
TAXON_NAME = 5
AUTHORS = 6
RANK = 7 ### Form, GENUS, InfraspecificName, SPECIES, Subform, SUBSPECIES, Subvariety, VARIETY
TAXONOMIC_STATUS = 8 ### Accepted, Artificial Hybrid, Homotypic_Synonym, Synonym, Unplaced
# ACCEPTED_KEW_ID = 9
ACCEPTED_NAME = 10
ACCEPTED_AUTHORS = 11
# PARENT_KEW_ID = 12
# PARENT_NAME = 13
# PARENT_AUTHORS = 14
# REVIEWED = 15
# PUBLICATION = 16
# ORIGINAL_NAME_ID = 17

### READ WCVP
ambiguous = {} ### name => True
synonymy = {} ### name => accepted name
with lzma.open(settings['wcvp'], mode = 'rt') as file:
	for k, line in enumerate(file):
		if k > 0:
			columns = line.rstrip().split('|')
			if columns[RANK] == 'SPECIES':
				if columns[TAXONOMIC_STATUS] in ('Accepted', 'Artificial Hybrid', 'Unplaced'):
					insertSpecies(ambiguous, columns[AUTHORS], columns[TAXON_NAME], columns[AUTHORS], columns[TAXON_NAME], synonymy)
				elif columns[TAXONOMIC_STATUS] in ('Homotypic_Synonym', 'Synonym', ''): ### '546041-1', '584168-1', '27083-2' have no status, but have accepted names (wtf?)
					insertSpecies(ambiguous, columns[AUTHORS], columns[TAXON_NAME], columns[ACCEPTED_AUTHORS], columns[ACCEPTED_NAME], synonymy)
				elif settings['debug'] == True:
					eprint(f"'{columns[KEW_ID]}' not stored! '{columns[TAXONOMIC_STATUS]}'")
			elif columns[RANK] in ('Form', 'Subform', 'SUBSPECIES', 'Subvariety', 'VARIETY'):
				if columns[TAXONOMIC_STATUS] in ('Accepted', 'Artificial Hybrid', 'Unplaced'):
					insertIntraspecific(ambiguous, columns[AUTHORS], columns[TAXON_NAME], columns[AUTHORS], columns[TAXON_NAME], synonymy)
				elif columns[TAXONOMIC_STATUS] in ('Homotypic_Synonym', 'Synonym'):
					insertIntraspecific(ambiguous, columns[AUTHORS], columns[TAXON_NAME], columns[ACCEPTED_AUTHORS], columns[ACCEPTED_NAME], synonymy)
				elif settings['debug'] == True:
					eprint(f"'{columns[KEW_ID]}' not stored!")
			elif columns[RANK] == 'InfraspecificName':
				if re.search(NOTHO, columns[TAXON_NAME]):
					if columns[TAXONOMIC_STATUS] in ('Accepted', 'Artificial Hybrid', 'Unplaced'):
						insertInfraspecificX(ambiguous, columns[AUTHORS], columns[TAXON_NAME], columns[AUTHORS], columns[TAXON_NAME], synonymy)
					elif columns[TAXONOMIC_STATUS] in ('Homotypic_Synonym', 'Synonym'):
						insertInfraspecificX(ambiguous, columns[AUTHORS], columns[TAXON_NAME], columns[ACCEPTED_AUTHORS], columns[ACCEPTED_NAME], synonymy)
					elif settings['debug'] == True:
						eprint(f"'{columns[KEW_ID]}' not stored!")
				else: ### convar., grex, lusus, microgene, modif., monstr., mut., nothosubsp., nothovar., proles, provar., stirps, sublusus, subproles, subso,
					if columns[TAXONOMIC_STATUS] in ('Accepted', 'Artificial Hybrid', 'Unplaced'):
						insertIntraspecific(ambiguous, columns[AUTHORS], columns[TAXON_NAME], columns[AUTHORS], columns[TAXON_NAME], synonymy)
					elif columns[TAXONOMIC_STATUS] in ('Homotypic_Synonym', 'Synonym'):
						insertIntraspecific(ambiguous, columns[AUTHORS], columns[TAXON_NAME], columns[ACCEPTED_AUTHORS], columns[ACCEPTED_NAME], synonymy)
					elif settings['debug'] == True:
						eprint(f"'{columns[KEW_ID]}' not stored!")
			elif columns[RANK] != 'GENUS' and settings['debug'] == True:
					eprint(f"'{columns[KEW_ID]}' not stored!")
if settings['debug'] == True:
	eprint(f"{len(synonymy):,} inputs stored, {len(ambiguous):,} ambiguous")



### REMOVE AMBIGUOUS WCVP
for key in ambiguous.keys():
	del synonymy[key]
if settings['debug'] == True:
	eprint(f"{len(synonymy):,} usable inputs stored")
	for key, value in synonymy.items():
		eprint(f"'{key}' => '{value}'")



### MAKE GENERIC KEY LISTS
genericSynonymy = {} ### genus => [keys]
for key in synonymy.keys():
	g = extractGenus(key)
	if g not in genericSynonymy:
		genericSynonymy[g] = []
	genericSynonymy[g].append(key)
if settings['debug'] == True:
	eprint(f"{len(genericSynonymy):,} genera stored")
	for key, value in genericSynonymy.items():
		eprint(f"{key} => {'|'.join(value)}")



### MATCH NAMES
with open(settings['input'], mode = 'rt', encoding = 'utf8', errors = 'replace') as file:
	for k, line in enumerate(file):
		columns = ftfy.fix_text(line).rstrip().split('\t')
		p = False
		if k == 0:
			insertColumns('algorithm', columns, f"corrected {columns[settings['column']]}")
			p = True
		else:
			if settings['debug'] == True:
				eprint(f"working on '{columns[settings['column']]}'...")
			genus = extractGenus(columns[settings['column']])
			tokens = tokenize(columns[settings['column']])
			if genus in genericSynonymy:
				for t, token in enumerate(tokens): ### exact token match
					if settings['debug'] == True:
						eprint(f"looking for exact matches to '{token}'...")
					if token in synonymy:
						insertColumns(f"exact token match ({int(100*((t+1)/len(tokens)))}%)", columns, synonymy[token])
						p = True
						if settings['debug'] == True:
							eprint(f"found '{token}' => '{synonymy[token]}'")
						break
			if p == False and settings['approximate'] != 'none': ### approximate token match
				names = []
				algorithm = 'approximate match within stated genus'
				if genus in genericSynonymy:
					names = genericSynonymy[genus]
					if settings['debug'] == True:
						eprint(f"trying agrep within {genus}...")
				else:
					count, score, matches = agrep(genericSynonymy.keys(), genus, 'all')
					if count >= 1:
						for match in matches:
							if match in genericSynonymy:
								names += genericSynonymy[match]
						algorithm = f"approximate match within approximated ({score}%) genus"
						if settings['debug'] == True:
							eprint(f"trying agrep within {', '.join(matches)}...")
				if len(names):
					for token in tokens:
						count, score, match = agrep(names, token, 'one')
						if count == 1:
							insertColumns(f"{algorithm} ({score}%)", columns, synonymy[match])
							p = True
							if settings['debug'] == True:
								eprint(f"agrep found '{match}' => '{synonymy[match]}'")
							break
		if p == True:
			print('\t'.join(columns))
		else:
			eprint(f"Could not resolve '{columns[settings['column']]}'!")