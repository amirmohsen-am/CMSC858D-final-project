import generate_query as gq

# refers to the number of clauses that will be in the queries
smallClauses = 5
mediumClauses = 20
largeClauses = 100

# these mins and maxs refer to how many genes will randomly be in a given clause
smallTermsMin = 5
smallTermsMax = 15

mediumTermsMin = 5
mediumTermsMax = 15

largeTermsMin = 15
largeTermsMax = 25

# total number of genes
numGenes = 9503

# number of queries that will go in each file
numQueries = 10

# E-CURD-21 is the name of the fly data files
f = open('E-CURD-21-random-queries-small-clauses-small-genes.in', 'w')
for i in range(numQueries):
	f.write(gq.randomQuery(numGenes, smallClauses, smallTermsMin, smallTermsMax) + "\n")
f.close()

f = open('E-CURD-21-random-queries-small-clauses-medium-genes.in', 'w')
for i in range(numQueries):
	f.write(gq.randomQuery(numGenes, smallClauses, mediumTermsMin, mediumTermsMax) + "\n")
f.close()

f = open('E-CURD-21-random-queries-small-clauses-large-genes.in', 'w')
for i in range(numQueries):
	f.write(gq.randomQuery(numGenes, smallClauses, largeTermsMin, largeTermsMax) + "\n")
f.close()


f = open('E-CURD-21-random-queries-medium-clauses-small-genes.in', 'w')
for i in range(numQueries):
	f.write(gq.randomQuery(numGenes, mediumClauses, smallTermsMin, smallTermsMax) + "\n")
f.close()

f = open('E-CURD-21-random-queries-medium-clauses-medium-genes.in', 'w')
for i in range(numQueries):
	f.write(gq.randomQuery(numGenes, mediumClauses, mediumTermsMin, mediumTermsMax) + "\n")
f.close()

f = open('E-CURD-21-random-queries-medium-clauses-large-genes.in', 'w')
for i in range(numQueries):
	f.write(gq.randomQuery(numGenes, mediumClauses, largeTermsMin, largeTermsMax) + "\n")
f.close()


f = open('E-CURD-21-random-queries-large-clauses-small-genes.in', 'w')
for i in range(numQueries):
	f.write(gq.randomQuery(numGenes, largeClauses, smallTermsMin, smallTermsMax) + "\n")
f.close()

f = open('E-CURD-21-random-queries-large-clauses-medium-genes.in', 'w')
for i in range(numQueries):
	f.write(gq.randomQuery(numGenes, largeClauses, mediumTermsMin, mediumTermsMax) + "\n")
f.close()

f = open('E-CURD-21-random-queries-large-clauses-large-genes.in', 'w')
for i in range(numQueries):
	f.write(gq.randomQuery(numGenes, largeClauses, largeTermsMin, largeTermsMax) + "\n")
f.close()

