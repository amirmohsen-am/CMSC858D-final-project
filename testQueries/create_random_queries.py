import generate_query as gq

def createQueryFiles():

	fewClauses = 10
	manyClauses = 100
	fewTermsMin = 5
	fewTermsMax = 15
	manyTermsMin = 90
	manyTermsMax = 110
	numGenes = 19000

	f = open('fewClauses_fewTerms', 'w')
	for i in range(5):
		f.write(gq.randomQuery(numGenes, fewClauses, fewTermsMin, fewTermsMax) + "\n")

	f.close()

	f = open('fewClauses_manyTerms', 'w')
	for i in range(5):
		f.write(gq.randomQuery(numGenes, fewClauses, manyTermsMin, manyTermsMax) + "\n")

	f.close()

	f = open('manyClauses_fewTerms', 'w')
	for i in range(5):
		f.write(gq.randomQuery(numGenes, manyClauses, fewTermsMin, fewTermsMax) + "\n")

	f.close()

	f = open('manyClauses_manyTerms', 'w')
	for i in range(5):
		f.write(gq.randomQuery(numGenes, manyClauses, manyTermsMin, manyTermsMax) + "\n")

	f.close()



