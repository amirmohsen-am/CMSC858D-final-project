from baselineCode import *
import sys, random
import json

# if len(sys.argv) < 2:
# 	print("Provide single-cell file")
# 	sys.exit()


class SingleCellSearch:
	# how many optimization we are using
	optimization_level = 3
	def __init__(self, fname, optimization_level=3, distance_file=None):
		self.optimization_level = optimization_level
		f = open(fname, 'r')
		f.readline()
		self.rows, self.cols, self.non_zero_cells = map(int, f.readline().split())
		self.numGenes = self.rows
		self.numCells = self.cols

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

		
		# Hamming Distance optimization
		if self.optimization_level >= 2:
			if distance_file is not None:
				self.load_distances(distance_file)
			else:
				print("Calculating hamming distance:")
				self.distance = [[0 for j in range(self.rows+1)] for i in range(self.rows+1)]
				for i in range(1, self.rows+1):
					if i%100 == 0:
						print("{:.2f}%".format(i/(self.rows+1)*100), end = " ")

					for j in range(i+1, self.rows+1):
						intersection = self.gene_cell[i].intersection(self.gene_cell[j])
						self.distance[i][j] = self.distance[j][i] = len(self.gene_cell[i]) + len(self.gene_cell[j]) - 2*len(intersection)
				print("hamming distance calculation done")



		print("Preprocessing done")
	
	def save_distances(self, distance_file):
		print("Saving distances to " + distance_file)
		with open(distance_file, 'w') as f:
			for i in range(self.rows+1):
				f.write(" ".join(map(str, self.distance[i])) + "\n")
		print("Saving done")
	
	def load_distances(self, distance_file):
		print("Loading distances from " + distance_file)
		self.distance = [[0 for j in range(self.rows+1)] for i in range(self.rows+1)]
		with open(distance_file, 'r') as f:
			i = 0
			for line in f.readlines():
				self.distance[i] = list(map(int, line.split(" ")))
				i += 1
		print("Loading done")


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

	# u and v could be negative
	def hamming_distance(self, u, v):
		if u*v > 0:
			return self.distance[u][v]
		else:
			return self.cols - self.distance[u][v]
			
		
	def search(self, query, optimization_level=3):
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
			intersected_cells = set(self.get_gene(smallest_gene))

			# look ahead optimization (removing results from representetive)
			if optimization_level >= 1:
				for x in list(intersected_cells):
					if x in result_cells:
						intersected_cells.remove(x)
			
			
			# hamming distance optimization
			if optimization_level >= 2:
				for i in range(2, len(priority)):
					if self.hamming_distance(priority[0][1], priority[1][1]) < self.hamming_distance(priority[0][1], priority[i][1]):
						priority[1], priority[i] = priority[i], priority[1]
				

			for i in range(1, len(priority)):
				v = priority[i][1]
				# make a copy of answer cells
				for x in list(intersected_cells):
					if v > 0:
						if not x in self.gene_cell[v]:
							intersected_cells.remove(x)
					else:
						if x in self.gene_cell[-v]:
							intersected_cells.remove(x)
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

def rand_neg(): 
	return random.randint(0, 1)*2-1

def create_meaninfulquery(scc_data, numClauses, numTerms):
	query = []
	numGenes = scc_data.numGenes
	for i in range(numClauses):
		clause = [random.randint(1, numGenes)]
		clause[0] *= rand_neg()

		closest = []
		for k in range(-numGenes, numGenes):
			if k != 0 and k != clause[0]:
				closest.append((scc_data.hamming_distance(clause[0], k), k))
		# print(len(closest))
		# print(len(closest[0]))
		closest.sort()
		for j in range(0, numTerms-1):
			clause.append(closest[j][1])
		query.append(clause)
	return query
		
def create_dumb_query(scc_data, numClauses, numTerms):
	query = []
	numGenes = scc_data.numGenes
	for i in range(numClauses):
		clause = [random.randint(-numGenes, numGenes)] * (numTerms-1)
		
		closest = []
		for k in range(-numGenes, numGenes):
			if k != 0 and k != clause[0]:
				closest.append((scc_data.hamming_distance(clause[0], k), k))
		closest.sort()

		clause.append(closest[-1][1])
		print(closest[-10:-1])
		# clause.append(random.randint(-numGenes, numGenes))
		
		query.append(clause)
	return query

			
