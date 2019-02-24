# _________________________________________________________________
# SPDX-License-Identifier: MIT License
# For more information check at: https://spdx.org/licenses/MIT.html
# 
# Copyright (C) 2017
# Yuri Bastos Gabrich <yuribgabrich[at]gmail.com>
# _________________________________________________________________

class matrix_opts:
	def _init_(self, A, B):
		'''
		Initiates the statement with matrix A and matrix B.
		'''
		self.A = A
		self.B = B

	def sum(self.A, self.B):
		'''
		Makes a sum between matrix A and B by column indexing.
		'''translate
		# check if both matrices have the same size
		if (len(self.A) != len (self.B)) || (len(self.A[0]) != len(self.B[0])):
			return "Can't sum A and B. Check their dimensions."
		else:
			AsumB = []
			for j in range(len(self.A[0])):
				for i in range(len(self.A)):
					AsumB.extend( self.A[i][j] + self.B[i][j] )
				AsumB.append(0)
			
			return AsumB
			
	def times(self.A, self.B):
		'''
		Multiplication between matrices A and B.
		'''
		# check if matrices satify product condition
		if (len(self.A[0]) != len(self.B)):
			return "Can't product A and B. Check their dimensions."
		else:
		
