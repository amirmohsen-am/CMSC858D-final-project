"""
To run the baseline, first format the matrix data by providing a .mtx file. This will return the number of genes, number of cells, 
and the matrix as a list of sets of cells (1 set per gene). Then you can call the lookUpQuery() or lookUpQueryTimed() function to 
get the results of the query. 

You can also look up bulk queries to get the how long it takes for all the queries to be looked up by calling the lookUpQueriesFromFile()
function.

This code is for doing baseline queries and timing them. If the queries are coming from a file then the format is 1 query per line,
with each query separating the clauses with ';' and the terms (genes) in those clauses separated by ','. NOT is represented by the 
gene being negative. For example a file with:

1,2,3;4,5,6
10,-11;12,13;14,15

would translate to the following queries:

(1 ^ 2 ^ 3) v (4 ^ 5 ^ 6)
(10 ^ NOT 11) v (12 ^ 13) v (14 ^ 15)

If using this file on the command line (without interactive python) then you can provide the path to the data file as the first command line argument and 
the path to the query file as the second command line argument. The total time it took to perform all queries in the query file will be output.

"""

import sys, time

# this function takes the name of at matrix file (.mtx) and returns the number of genes, number of cells, and the matrix as a list of sets (1 set of cells per gene) 
def formatData(filename):

	matrix = []

	f = open(filename, "r")
	f.readline()
	stats = f.readline().split(" ")
	numGenes = int(stats[0])
	numCells = int(stats[1])
	# the 0 row will be empty, but it makes it less confusing because the data starts at 1, not 0
	for i in range(int(numGenes) + 1):
		matrix.append(set())

	for line in f.readlines():
		if len(line) > 0:
			vals = line.split(" ")
			matrix[int(vals[0])].add(int(vals[1]))

	f.close()

	return numGenes, numCells, matrix


# the query will be a list of lists, where the inner lists represent "AND"s and the inner lists are "OR"ed together (DNF)
def lookUpQuery(numGenes, numCells, matrix, query):

	results = set()
	# there is no cell 0
	for i in range(1, numCells + 1):
		for disjunction in query:
			disjunctionVerdict = True
			for j in disjunction:
				if j < 0 and i in matrix[-j]:
					disjunctionVerdict = False
					break
				elif j > 0 and i not in matrix[j]:
					disjunctionVerdict = False
					break
			# if all parts of disjunction were checked and none were false, add this cell to the results
			if disjunctionVerdict:
				results.add(i)

	return results


# this function calls the times how long it took to look up a single query
def lookUpQueryTimed(numGenes, numCells, matrix, query):

	timeBefore = time.perf_counter()
	results = lookUpQuery(numGenes, numCells, matrix, query)
	timeAfter = time.perf_counter()

	return results, timeAfter - timeBefore


# this function is for looking up queries that are coming from a file and timing them (the file should just have one query per line)
def lookUpQueriesFromFile(numGenes, numCells, matrix, queryFile):

	totalTime = 0.0

	queries = getQueriesFromFile(queryFile)
	for query in queries:
		results, queryTime = lookUpQueryTimed(numGenes, numCells, matrix, query)

		totalTime += queryTime

	return totalTime


# this function returns a list of queries from a file (one query per line in the file. The clauses are separated by ';' and the terms within the clauses are separated by ','.)
def getQueriesFromFile(filename):

	queries = []

	f = open(filename, "r")
	for line in f.readlines():
		query = []

		# split the line into individual clauses
		clause_strings = line.split(";")
		for clause_string in clause_strings:
			clause = []

			# split the clause into the individual terms (genes)
			terms = clause_string.split(",")
			for term in terms:
				clause.append(int(term))

			query.append(clause)

		queries.append(query)

	f.close()

	return queries


# this function is just for testing my lookUpQuery() function
def runTests(queriesFile, dataFile):

	queries = getQueriesFromFile(queriesFile)
	print("Loading matrix data...")
	numGenes, numCells, matrix = formatData(dataFile)

	count = 0
	for query in queries:
		print("Current query:",query)
		expectedResults = getExpectedResults(numCells, matrix, query)
		baselineResults = lookUpQuery(numGenes, numCells, matrix, query)

		if expectedResults != baselineResults:
			print("Query Number:",count," FAILED")
		else:
			print("Query Number:",count," passed")

		count += 1


# this function is for testing the code. It basically just gets the results by directly looking at the matrix. 
# The returned value can then be compared to the results from the lookUpQuery() function
def getExpectedResults(numCells, matrix, query):

	results = set()

	for clause in query:
		clauseResults = set()
		first = True
		for gene in clause:
			if first:	
				clauseResults.update(getGene(numCells, matrix, gene))
				first = False
			else:
				clauseResults = clauseResults.intersection(getGene(numCells, matrix, gene))

		results.update(clauseResults)

	return results

# this function returns either the cells that are in the set or the cells NOT in the set depending on the parity of the gene
def getGene(numCells, matrix, gene):

	if gene > 0:
		return matrix[gene]
	else:
		notGene = []
		for i in range(1, numCells + 1):
			if i not in matrix[-gene]:
				notGene.append(i)
		return notGene


if __name__ == "__main__":

	dataFile = sys.argv[1]
	queryFile = sys.argv[2]

	print("Formatting data...")
	numGenes, numCells, matrix = formatData(dataFile)
	
	print("Doing the search...")
	totalTime = lookUpQueriesFromFile(numGenes, numCells, matrix, queryFile)

	print("Total time was: ",totalTime)