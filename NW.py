#!/usr/bin/env python3
"""
Author: Andreas Chrysostomou
Code to implement the Needleman – Wunsch algorithm.
"""

#functions
def parse_sequences(filename):
	"""
	The function parses the sequences from the fasta file.
	:param filename: string, the name of the file containing the sequences
	:return sequences: tuple of strings
	"""
	data = []
	f = open(filename, 'r')
	for line in f:
		if not line.strip():
			continue
		if line.startswith('>'):
			continue
		else:
			sequence = line.strip()
			data.append(sequence)
	f.close()
	return data

def parse_score_matrix(filename):
	"""
	The function parses the score matrix from the file.
	:param filename: string, the name of the file containing the score matrix
	:return score_matrix: tuple, a tuple in form (reference, matrix) where score is a dictionary and reference a list
	to guide in which position to look for.
	"""
	f = open(filename, 'r')
	reference = f.readline().strip().split()
	matrix = {}
	for line in f:
		line = line.strip().split()
		key = line[0] # The first element of the line is the amino acid or nucleotide that we will use as a key in our
		# dictionary and the rest ([1:]) are the values that correspond the elements of the reference list.
		# We use the map function to convert the elements of the list into integers all at once.
		scores = list(map(int, line[1:]))
		matrix[key] = scores
	f.close()
	score_matrix = (reference, matrix)
	return score_matrix

def get_score(elem1, elem2):
	"""
	The function returns the score for the corresponding elements.
	:param elem1: string, a nucleotide or amino acid (element) from the first sequence
	:param elem2: string, a nucleotide or amino acid (element) from the second sequence
	:return score: int, the corresponding score of the elements
	"""
	reference, matrix = SCORING_SYSTEM_MATRIX
	# This should be used as global variable to avoid the parsing or passing the score matrix
	# to the lookup function every time we need to do a calculation.
	index = reference.index(elem1)
	key = elem2
	score = matrix[key][index]
	return score


def initiate_matrices(seq1, seq2, gap):
	"""
	The function initiates the matrices for the calculations.
	:param seq1: string, the first sequence to be aligned
	:param seq2: string, the second sequence to be aligned
	:param gap: gap, int, the gap penalty
	:return matrices: tuple, (calc_matrix, symb_matrix) contains the matrices (essentially lists of lists)
	for the following calculations
	"""
	rows = len(seq1) + 1 # +1 is added for the gap row.
	cols = len(seq2) + 1 # +1 is added for the gap column.

	# calc_matrix will carry the calculations of the scores
	calc_matrix = [[0 for i in range(cols)] for j in range(rows)]
	# Fill in the gap row and column.
	for j in range(1,cols):
		calc_matrix[0][j] = j * gap
	for i in range(1,rows):
		calc_matrix[i][0] = i * gap

	# symb_matrix will carry the symbols for the trace back
	# '<' means an indel. Move leftwards.
	# '^' means an indel. Move upwards.
	# '*' means a match or mismatch. MOve diagonally.
	symb_matrix = [['0' for i in range(cols)] for j in range(rows)]
	# Fill in the gap row and column.
	for j in range(1,cols):
		symb_matrix[0][j] = '<'
	for i in range(1,rows):
		symb_matrix[i][0] = '^'

	# Pack the matrices into a tuple to return.
	matrices = (calc_matrix, symb_matrix)
	return matrices

def fillup_matrices(seq1, seq2, gap):
	"""
	The function will full up all the cells in the initiated matrices. It runs the actual algorithm.
	:param seq1: string, the first sequence to be aligned
	:param seq2: string, the second sequence to be aligned
	:param gap: gap, int, the gap penalty
	:return matrices: tuple, (calc_matrix, symb_matrix) contains the matrices (essentially lists of lists)
	"""
	# Call the initiate_matrices function to prepare the matrices
	initiated_matrices = initiate_matrices(seq1=seq1, seq2=seq2, gap=gap)
	calc_matrix, symb_matrix = initiated_matrices
	# We add a whitespace to align the indexes of the sequence strings with the matrix coordinates. This is needed
	# to take care for the gap column and line.
	sequence1 = ' ' + seq1
	sequence2 = ' ' + seq2
	# We fill up all the cells.
	for i in range(1,len(sequence1)):
		for j in range(1,len(sequence2)):
			score = get_score(sequence1[i], sequence2[j])
			diag = calc_matrix[i-1][j-1] + score
			diag = (diag,'*')
			up = calc_matrix[i-1][j] + gap
			up = (up,'^')
			left = calc_matrix[i][j - 1] + gap
			left = (left, '<')
			maximum = max([diag, up, left])
			max_score, max_symbol = maximum
			calc_matrix[i][j] = max_score
			symb_matrix[i][j] = max_symbol
	# Pack the matrices into a tuple to return.
	matrices = (calc_matrix, symb_matrix)
	return matrices


def construct_alignment(seq1, seq2, symb_matrix):
	"""
	The function creates the alignment based on the symbol matrix we created before.
	:param seq1: string, the first sequence to be aligned
	:param seq2: string, the second sequence to be aligned
	:param symb_matrix: list of lists, the matrix containing the symbols needed for the traceback and alignment
	construction
	:return final_alignment: tuple of strings, contains the final alignments
	"""
	# We add a whitespace to align the indexes of the sequence strings with the matrix coordinates. This is needed
	# to take care for the gap column and line.
	sequence1 = ' ' + seq1
	sequence2 = ' ' + seq2
	alignment1 = ''
	alignment2 = ''
	i = len(seq1)
	j = len(seq2)
	while i > 0 and j >0:
		if symb_matrix[i][j] == '*':
			alignment1 = sequence1[i] + alignment1
			alignment2 = sequence2[j] + alignment2
			i = i - 1
			j = j - 1
		elif symb_matrix[i][j] == '^':
			alignment1 = sequence1[i] + alignment1
			alignment2 = '-' + alignment2
			i = i-1
		elif symb_matrix[i][j] == '<':
			alignment1 = '-' + alignment1
			alignment2 = sequence2[j] + alignment2
			j = j - 1

	#In the case we reach i=0 or j=0
	while i > 0:
		alignment1 = sequence1[i] + alignment1
		alignment2 = '-' + alignment2
		i = i - 1

	while j > 0:
		alignment1 = '-' + alignment1
		alignment2 = sequence2[j] + alignment2
		j = j - 1

	final_alignment = (alignment1, alignment2)
	return final_alignment

def output(final_alignment, matrices):
	"""
	The function prints the output.
	:param final_alignment: tuple of strings, contains the final alignments
	:param matrices: tuple, (calc_matrix, symb_matrix) contains the matrices (essentially lists of lists)
	:return None
	"""
	calc_matrix, symb_matrix = matrices
	print('The final alignment is:')
	for alignment in final_alignment:
		print(alignment)

	print('Here is the matrix with calculations:')
	for row in calc_matrix:
		print(row)
	print('Here is the traceback matrix:')
	for row in symb_matrix:
		print(row)
	return None

if __name__ == '__main__':
	FASTA_FILENAME = 'Fasta.txt'
	DATA = parse_sequences(filename= FASTA_FILENAME)
	SEQ1, SEQ2 = DATA
	SCORING_SYSTEM_FILENAME = 'BLOSUM62.txt'
	SCORING_SYSTEM_MATRIX = parse_score_matrix(filename = SCORING_SYSTEM_FILENAME)
	GAP = -1
	MATRICES = fillup_matrices(seq1=SEQ1, seq2=SEQ2, gap=GAP)
	CALC_MATRIX, SYMB_MATRIX = MATRICES
	FINAL_ALIGNMENT = construct_alignment(seq1=SEQ1, seq2=SEQ2, symb_matrix= SYMB_MATRIX)
	output(final_alignment=FINAL_ALIGNMENT, matrices=MATRICES)
  
