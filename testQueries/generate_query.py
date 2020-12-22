import random
def randomQuery(numGenes, numClauses, minClauseSize, maxClauseSize):
    query = ""
    for clause in range(numClauses):
        clauseSize = random.randint(minClauseSize, maxClauseSize)
        for atom in range(clauseSize):
            gene = random.randint(-numGenes, numGenes)
            while gene == 0:
                gene = random.randint(-numGenes, numGenes)
            query += str(gene) + ","
        query = query[:-1] + ";"
    return query
