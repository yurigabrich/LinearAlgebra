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
        aug_matrix = []
        for i in range(len(lines)):
            aug_matrix.append(lines[i].split()) #column delimiter ' ' --> single space
                
            # converting a row of strings to a row of floats
            for elem in range(len(aug_matrix[i])):
                aug_matrix[i][elem] = float(aug_matrix[i][elem])
                
        # column: size of the first row
        cols = len(aug_matrix[0])

        # cheking if some element was forgotten
        for k in range(len(aug_matrix)):
            if len(aug_matrix[k]) != cols:
                raise IndexError
        
        # end of load_matrix, if everything is right
        print('\n', 'Matrix', len(aug_matrix), 'x', cols, 'loaded.')
        in_file.close()
        return aug_matrix, cols
        
    except IndexError:
        missed_element = {} # dictionary{ cols: rows }
        for n in range(len(aug_matrix)):
            missed_element[len(aug_matrix[n])] = n
        
        issued_row = missed_element.get(min(missed_element.keys()))
        
        print('\n',"Ops! Apparently you've eaten some element(s) of the matrix.")
        print(" Check row", str(1 + issued_row)+":")
        print('    ', aug_matrix[issued_row])
        
    except ValueError:
        print('\n',"Hey, only numbers! Check row", str(i+1))
        
    #except SyntaxError: #seria parenteses, colchetes...


def ncr(n, k): #talvez seja apagado
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
    
    matrix (list of list of floats): matrix organized by rows
    
    Returns:
            - boolean indicating initial condition about square behaviour
            - list of list with all square submatrices with the highest order
    '''
    
    # checking if it is a square matrix
    if len(matrix) == len(matrix[0]):
        return True, matrix
    else:
        # looking for the number of submatrices we can get
        min_index = min(len(matrix),len(matrix[0]))
        
        submatrices = itertools.combinations(np.transpose(matrix), min_index)
        
        return False, submatrices


def det_matrix(check, submatrices):
    '''
    Calculates the determinants of submatrices accordingly with the
    previous information about dimension of the main matrix (square or not).
    
    Returns (int): the order of the submatrices
    '''
    
    if check:
        det = np.linalg.det(submatrices)
        order = len(list(submatrices))
    else:
        det = 0
        for k in submatrices:
            det += np.linalg.det(k)
            if k != 0:
                order = len(list(submatrices))
                break
    
    return order


def solutionize(matrix, i, j):
    '''
    Indicates one of the three types of a solution in linear algebra system.
    
    matrix (list of list of floats): matrix organized by rows
    i (int): number of rows of the matrix 
    j (int): number of columns of the matrix 
    
    Returns:
            - boolean statement to indicate if the echelon must be made,
            - string with type of solution
    '''
    # first, find A
    A = []
    for k in range(i):
        A.append(matrix[k][:-1])
        
    # calculates the determinante of A and get the rank(A)
    check_A, submatrix_A = is_square(A)
    p_A = det_matrix(check_A, submatrix_A)
    
    # calculates the determinante of M and get the rank(M)
    check_M, submatrix_M = is_square(matrix)
    p_M = det_matrix(check_M, submatrix_M)
    
    # compare rank of A and M just once
    postos = (p_A == p_M)
    
    # identify solution type
    if p_A == min(i,j-1): # j-1 = number of columns of matrix A
        if i == j-1:
            return True, "SPD"
        elif i < j-1:
            return True, "SPI"
        else: # i > j-1
            if postos:
                return True, "SPD"
            else:
                return False, "SI"
    else: #p_A < min{i,j-1}
        if postos:
            return True, "SPI"
        else:
            return False, "SI"


def echelon(aug_matrix, i, j, stop = False):
    '''
    Does the echelon of a matrix of n-size, recursively.
    Put the greater element as pivot for each row.
    Start to do echelon from top to bottom.
    
    Inputs:
        - aug_matrix (list of list of floats): full matrix organized by lines
        - i (int): number of rows of aug_matrix
        - j (int): number of columns of aug_matrix
    
    Returns: stepped matrix organized as a list of list of floats
    '''
    if not stop:
        # identify pivot element and position
        i_turn = len(aug_matrix) - i
        j_turn = len(aug_matrix[0]) - j
        pivot = aug_matrix[i_turn][j_turn]
        
        # arrange rows positions accordingly with non-null initial pivot
        for n in range(i_turn + 1, i):
            #change line positions if needed        
            if aug_matrix[n][j_turn] > pivot:
                temp_row = aug_matrix[n].copy()
                aug_matrix[n] = aug_matrix[i_turn].copy()
                aug_matrix[i_turn] = temp_row.copy()
                # UPDATE PIVOT
                pivot = aug_matrix[i_turn][j_turn]
        
        # look for values non-null on the left of pivot item and clean them
        for k in range(j_turn):
            if (aug_matrix[i_turn][k]) != 0.00:
                # change the value of the entirely row
                num = aug_matrix[i_turn][k]
                denom = aug_matrix[k][k]
                for x in range(len(aug_matrix[0])): # because 'j' changes recursively!
                    aug_matrix[i_turn][x] -= abs( aug_matrix[k][x] * (num / denom))
            # UPDATE PIVOT
            pivot = aug_matrix[i_turn][j_turn]
        
        # make pivot element be a unit,
        # as well as change the value of the entirely row
        for k in range(len(aug_matrix[0])): # because 'j' changes recursively!
            aug_matrix[i_turn][k] = round(aug_matrix[i_turn][k] / pivot, 2)
        
        #stop criterion
        if (min(i_turn,j_turn)-1) == 0:
            stop = True
        
        return echelon(aug_matrix, i-1, j-1, stop)
    
    else:
        return aug_matrix #stepped


def calc_var(row, pos, variables_x):
    '''
    Finds the value of the unkown by simple math accordingly with row specs
    
    - row (list of floats): desired line of the echelon matrix
    - pos (int): column index of the unkown
    - variables_x (list of floats): value of the unknowns
    
    Returns: value of the unkown (float)
    '''
    value = row[-1]

    #print('\n', value, variables_x, end='')

    for k in range(len(row)-1, pos, -1):
        
        value -= (row[k] * variables_x[len(row) -1 -k])

        #print('', value, k, '\n')

    return round(value, 2)


#-------------------------------------------------------------------
class Gauss(object):
    
    def __init__(self, matrix):
        '''
        Initializes Gauss object. The only value needed is the matrix.
        
        a Gauss object has four attributes:
            self.matrix (list of string, the file selected on MATRIX)
            self.aug_matrix (list of list of floats, the matrix above converted to floats and organized by rows)
            self.cols (total number of columns in the matrix)
            self.memory_aug_matrix (a copy of inputed matrix, because of memory address)
            
        Returns NOTHING!
        '''
        self.matrix = matrix
        # read input file with unknown dimension of an augmented matrix (A)
        self.aug_matrix, self.cols = load_matrix(self.matrix)
        self.memory_aug_matrix = copy.deepcopy(self.aug_matrix)
        return None
    
    
    def get_initial_matrix(self):
        '''
        Used to safely access self.matrix outside of the class
        
        Returns: NOTHING! Only print the original matrix.
        '''
        for x in range(len(self.memory_aug_matrix)):
            print(self.memory_aug_matrix[x])
        return None
        
        
    def get_size(self):
        '''
        Returns: the number of rows (i) and columns (j) of the aug_matrix
        '''
        return len(self.aug_matrix), self.cols
    
        
    def result(self):
        '''
        Makes the Gauss echelon method
        
        Returns: value of the variables (por enquanto matriz escalonada somente)
                DEFINIR O QUE RETORNAR AINDA!
        '''
        
        # gets the real size of the augmented matrix, not index of the 'list' argument
        i, j = self.get_size()
        
        # verify solution
        do_gauss, sys_type = solutionize(self.aug_matrix, i, j)
        
        # do echelon
        if do_gauss:
            step_matrix = echelon(self.aug_matrix, i, j)
            variables_x = [0]*(j-1) #(j-1) = number of columns of matrix A
            
            if sys_type == "SPD":
                for k in range(i, 0, -1):
                    variables_x[k-1] = calc_var(step_matrix[k-1], k-1, variables_x)
                    #print(variables_x)
                
            else: #sys_type == "SPI"            
                # do variables substitution (recursive form)
                    # if (eq < var) --> print lambda
                        # (i , 0, -1)
                    # else
                #calc_var(step_matrix[0],0)
                # update set of variables
                for k in range(i):                    
                    variables_x.append(step_matrix[k-1][k-1])    
            
            # printing final solution
            print("", sys_type, "=", variables_x)
            print("echelon form:")
            for x in range(len(step_matrix)):
                print(step_matrix[x])
        
        else:
            print("", sys_type)
        
        return None

#-------------------------------------------------------------------
M = 'matrix_SPD.txt'
#M = 'matrix_SPI.txt'
#M = 'matrix_SI.txt'
#M = 'matrixB.txt'

test = Gauss(M)
test.get_initial_matrix()
print('\n', "Gauss result:", end="")
test.result()
#print(test.get_size())
