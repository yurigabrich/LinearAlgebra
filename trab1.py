'''
Created on Wed Jan 18 11:30:52 2017

Trabalho 1 - Álgebra Linear
Prof. Zochil Arenas
Auno: Yuri Bastos Gabrich
'''

#import string
#import sys
import copy
import numpy as np
import functools
import operator as op
import itertools

def load_matrix(file_name):
    '''
    file_name (string): the name of the file containing 
    the list of string elements to load    
    
    Returns: a list of list of valid float elements, that is a list of rows
             and the number or columns, that'll be used to stop iterations
    
    Depending on the size of the matrix, this function may
    take a while to finish.
    ''' 
    print('\n''Loading matrix from file...')
    # inFile: file
    in_file = open(file_name, 'r')
    
    # lines: STRING
    lines = in_file.readlines() # line delimiter \n
    
    try:
        rows = []
        for i in range(len(lines)):
            rows.append(lines[i].split()) #column delimiter ' ' --> single space
                
            # converting a row of strings to a row of floats
            for elem in range(len(rows[i])):
                rows[i][elem] = float(rows[i][elem])
                
        # column: size of the first row
        cols = len(rows[0])

        # cheking if some element was forgotten
        for k in range(len(rows)):
            if len(rows[k]) != cols:
                raise IndexError
        
        # end of load_matrix, if everything is right
        print('\n', 'Matrix', len(rows), 'x', cols, 'loaded.')
        in_file.close()
        return rows, cols
        
    except IndexError:
        missed_element = {} # dictionary{ cols: rows }
        for n in range(len(rows)):
            missed_element[len(rows[n])] = n
        
        issued_row = missed_element.get(min(missed_element.keys()))
        
        print('\n',"Ops! Apparently you've eaten some element(s) of the matrix.")
        print(" Check row", str(1 + issued_row)+":")
        print('    ', rows[issued_row])
        
    except ValueError:
        print('\n',"Hey, only numbers! Check row", str(i+1))
        
    #except SyntaxError: #seria parenteses, colchetes...


def ncr(n, k):
    '''
    The following program calculates nCr in an efficient manner
    (compared to calculating factorials etc.)
    fonte: http://stackoverflow.com/a/4941932/6770397
    
    Inputs:
        - n = number of elements
        - k = k-combinations of n
    
    Returns (int): Combination of n and k
    '''
    r = min(k, n-k)
    if r == 0:
        return 1
    else:
        numer = functools.reduce(op.mul, range(n, n-r, -1))
        denom = functools.reduce(op.mul, range(1, r+1))
        return numer//denom


def is_square(matrix):
    '''
    Check if the matrix is square
    
    matrix (list of list of floats): matrix organized by lines (the 'rows' loaded)
    
    Returns all square submatrices with the highest order
    '''
    
    # checking if it is a square matrix
    if len(matrix) == len(matrix[0]):
        return matrix
    else:
        #submatrices = []

        # looking for the number of submatrices we can get
        min_index = min(len(matrix),len(matrix[0]))
        max_index = max(len(matrix),len(matrix[0]))
        arrange = ncr(max_index, min_index) #talvez não precise mais disso!!
        print(arrange)
        c = itertools.combinations(np.transpose(matrix), 3)
        count = 1
        for k in c:
            print(count)
            print(np.transpose(k))
            count+=1
        
        #for x in range(arrange):
        #    submatrices.append()
            
        return arrange, c


def solutionize(matrix, i, j):
    '''
    Indicates one of the three types of a solution in linear algebra system.
    
    matrix (list of list of floats): matrix organized by lines (the 'rows' loaded)
    
    Returns:
            - boolean statement to indicate if the echelon must be made,
            - string with type of solution
    '''
    A = []
    
    # checking if A is a square matrix
    if len(A) == len(A[0]):
        p_A = np.linalg.det(A)
    else:
        for x in range(len(matrix)):
            A.append(matrix[x][:-1]) # eliminating the last column
        #criar uma função para diminuir as colunas!!!
    
    # identify solution type
    if p_A == i:
        return True, "SPI"
    else:
        #p_M = numpy.linalg.det(matrix)
        
        # compare rank of A and M just once
        postos = False#(p_A != p_M)
        
        if p_A == j:
            if postos:
                return False, "SI"
            else:
                return True, "SPD"
        else: # must be only p_A < min{i, j}
            if postos:
                return False, "SI"
            else:
                return True, "SPI"


def echelon(rows, cols, i, j, stop = False):
    '''
    Does the echelon of a matrix of n-size, recursively
    
    rows (list of list of floats): full matrix organized by lines
    cols (int): number of columns on matrix
    i (int): number of rows
    j (int): number of columns
    
    Returns: stepped matrix organized as a list of list of floats
    '''
    if not stop:
        # identify pivot (PIV)
        pivj = cols - j
        pivi = len(rows) - i
        PIV = rows[pivi][pivj]
        
        # look for the greatest (little) item on the pivot column, and define the
        # corresponding row as main row (maybe because of computacional error problem)   
        for n in range(pivi + 1, len(rows)):
            #change line positions if needed        
            if rows[n][pivj] > PIV:
                temp_row = rows[n].copy()
                rows[n] = rows[pivi]          
                rows[pivi] = temp_row
                PIV = rows[pivi][pivj] #update pivot
                
        # look for values non-null on the left of pivot item and clean then
        for x in range(pivj):
            if (rows[pivi][pivj - x]) != 0.00:
                # change the value of the entirely row
                X = rows[pivi][x]
                for n in range(pivj):
                    rows[pivi][n] -= X * abs(rows[x][n])
            PIV = rows[pivi][pivj] #update pivot

        # make pivot element be a unit,
        # as well as change the value of the entirely row
        for n in range(cols):
            rows[pivi][n] = round(rows[pivi][n] / PIV, 2)
        
        #stop criterion
        if (min(i,j)-1) == 0:
            stop = True
        
        return echelon(rows, cols, i-1, j-1, stop)
    
    else:
        return rows #echeloned


#-------------------------------------------------------------------
class Gauss(object):
    
    def __init__(self, matrix):
        '''
        Initializes Gauss object. The only value needed is the matrix.
        
        a Gauss object has four attributes:
            self.matrix (list of string, the file selected on MATRIX)
            self.rows (list of list of floats, the matrix above converted to floats and organized by rows)
            self.cols (total number of columns in the matrix)
            self.memory_rows (a copy of inputed matrix, because of memory address)
            
        Returns NOTHING!
        '''
        self.matrix = matrix
        # read input file with unknown dimension of an augmented matrix (A)
        self.rows, self.cols = load_matrix(self.matrix)
        self.memory_rows = copy.deepcopy(self.rows)
        return None
    
    
    def get_initial_matrix(self):
        '''
        Used to safely access self.matrix outside of the class
        
        Returns: NOTHING! Only print the original matrix.
        '''
        for x in range(len(self.memory_rows)):
            print(self.memory_rows[x])
        return None
        
        
    def get_size(self):
        '''
        Returns: the number of rows (i) and columns (j) of the matrix
        '''
        is_square(self.rows)
        
        return #len(self.rows), self.cols
        
    
    # how to identify the 'posto(A)'?
    
        
    def result(self):
        '''
        Makes the Gauss echelon method
        
        Returns: value of the variables (por enquanto matriz escalonada somente)
                DEFINIR O QUE RETORNAR AINDA!
        '''
        
        # do echelon (recursive form)
        i, j = self.get_size()
        
        # verify solution
        go, sys_type = solutionize(self.rows, i, j)
        
        if go:
            stepped_matrix = echelon(self.rows,self.cols,i,j)
            
            #print()
            #for x in range(len(stepped_matrix)):
            #    print(stepped_matrix[x])
            
            # do variables substitution (recursive form)
            
                # if (eq < var) --> print lambda
                    
                # else
            
            for x in range(len(stepped_matrix)):
                print(stepped_matrix[x])
        else:
            print(sys_type)
        
        return None

#-------------------------------------------------------------------
MATRIX = 'matrixA.txt'
#MATRIX = 'matrix-A.txt'
#MATRIX = 'matrixB.txt'
#MATRIX = 'matrixC.txt'
test = Gauss(MATRIX)
#print(test.get_initial_matrix())
#print('\n', "Gauss result:")
#print(test.result())
print(test.get_size())
