from baselineCode import *
import sys

# if len(sys.argv) < 2:
# 	print("Provide single-cell file")
# 	sys.exit()


class SingleCellSearch:
	def __init__(self, fname):
		f = open(fname, 'r')
		f.readline()
		self.rows, self.cols, self.non_zero_cells = map(int, f.readline().split())

		# genes are 1 based
		
		self.cnt = [0 for i in range(self.rows+1)]

		# initalize empty set of cells for each gene
		self.gene_cell = [set() for i in range(self.rows+1)]

		self.small = 0

		print("Preprocessing started:")
		for q in range(self.non_zero_cells):
			if q % 500000 == 0:
				print("{:.2f}%".format(q/self.non_zero_cells*100), end = " ")
			s = f.readline().split()
			i, j, c = int(s[0]), int(s[1]), float(s[2])
			#


			if c < 0.1:
				self.small += 1
	#m[i].append(j);
			self.gene_cell[i].add(j)
			self.cnt[i] = self.cnt[i]+1

		print(self.non_zero_cells)
		print(float(self.non_zero_cells)/(self.rows*self.cols))


		self.srt = [(self.cnt[i], i) for i in range(self.rows)]

		self.srt.sort(reverse=True)

		print("Preprocessing done")

	# v < 0: not of v
	def get_gene(self, v):
		assert(v != 0)
		if v > 0:
			return self.gene_cell[v]
		if v < 0:
			not_g = []
			for i in range(1, self.cols+1):
				if i not in self.gene_cell[-v]:
					not_g.append(i)
			return not_g
		
	def search(self, query):
		result_cells = set()
		for clause in query:
			priority = []
			for v in clause:
				n_cells = 0
				if (v < 0):
					n_cells = self.rows - len(self.gene_cell[-v])
				else:
					n_cells = len(self.gene_cell[v])
				priority.append((n_cells, v))

			priority.sort()
			
			smallest_gene = priority[0][1]
			intersected_cells = self.get_gene(smallest_gene)
			
			for i in range(1, len(priority)):
				v = priority[i][1]
				# make a copy of answer cells
				tmp = set(intersected_cells)
				for x in intersected_cells:
					if v > 0:
						if not x in self.gene_cell[v]:
							tmp.remove(x)
					else:
						if x in self.gene_cell[-v]:
							tmp.remove(x)
				intersected_cells = tmp
			# query_str = '^'.join("g_"+str(v) for v in clause)			
			# print("clause:", query_str)
			# print(len(intersected_cells))
			result_cells.update(intersected_cells)
		return result_cells
		



# this function is just for testing my lookUpQuery() function
def test_code(scc_data, baseline_data, queries):

	numGenes, numCells, matrix = baseline_data

	count = 0
	for query in queries:
		#print("Current query:",query)
		baselineResults = lookUpQuery(numGenes, numCells, matrix, query)
		code_results = scc_data.search(query)


		if code_results != baselineResults:
			print("Query Number:",count," FAILED")
			print("Difference:", code_results.symmetric_difference(baselineResults))
		else:
			print("Query Number:",count," passed")

		count += 1