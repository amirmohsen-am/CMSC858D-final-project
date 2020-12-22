from SingleCellSearch import SingleCellSearch

#fname = "/Users/mos/Desktop/bio-project/ebi.ac.uk/E-MTAB-7051-quantification-raw-files/E-MTAB-7051.aggregated_filtered_counts.mtx"
fname = "../human_dataset/E-MTAB-9221.aggregated_filtered_counts.mtx"


data = SingleCellSearch(fname)

n_query = int(input("Number of queries:"))
# genome input is given 1 based from ebi.ac.uk

# queries are in a form of disjunctive normal form (Or of And's)
# first line number of clauses are given
# clauses are given in each line, on each line the genes are given (genes to be 'and')
for q in range(n_query):
	n_clause = int(input("Number of clauses:"))
	query = []
	for c in range(n_clause):
		# genes to be intersected
		clause = list(map(int, input().split()))
		query.append(clause)
	result = data.search(query)
	print("cells found: {}".format(len(result)))




