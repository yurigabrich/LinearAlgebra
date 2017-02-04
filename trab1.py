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
#import functools
#import operator as op
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


def split_matrix(matrix, number):
    '''
    Divides a matrix in several submatrices accordingly to defined number.
    Returns: a list of submatrices
    '''
    return itertools.combinations(np.transpose(matrix), number)


def squares(matrix):
    '''
    Check if the matrix is square and returns a list of list with all
    square submatrices with the highest order.
    
    matrix (list of list of floats): matrix organized by rows
    '''
    # checking if it is a square matrix
    if len(matrix) == len(matrix[0]):
        return [matrix]
    else:
        # looking for the number of submatrices we can get
        min_index = min(len(matrix),len(matrix[0]))
        return split_matrix(matrix, min_index)


def catch_zero(submatrix, dim):
    '''
    Looks for null, equal or proportional rows (or columns).
    
    - submatrix = biggest square submatrix
    - dim = dimension of the square submatrix

    Returns: boolean statement indicating if the determinant will be zero.
    '''
    for k in range(dim):
            # look for null rows
            count_zero = submatrix[k].count(0)
            if count_zero == dim:
                # det = 0 DETECTED!
                return True
                
            for s in range(k+1, dim, 1):
                # look for equal rows
                if submatrix[k] == submatrix[s]:
                    # det = 0 DETECTED!
                    return True
            
                # look for ratio number between rows
                ratio = submatrix[k][0] % submatrix[s][0]
                remainder = 0
                for t in range(k+1, dim, 1):
                    # if one is multiple of another
                    if (submatrix[s][t] == ratio * submatrix[k][t]):
                        remainder += 1
                if remainder == dim:
                    # det = 0 DETECTED!
                    return True
    
    return False # det != 0


def indirect_det(submatrix):
    '''
    Indicates if the determinant of a submatrix (can be a list of submatrices)
    will be zero or not, and returns the rank of the submatrix.
    '''
        
    dim = len(submatrix[0])
    if dim != 1:
        count_zero_det = 0
        
        for square_submatrix in submatrix:
            Tsubmatrix = np.transpose(square_submatrix)
            
            is_det_zero = catch_zero(square_submatrix, dim) or catch_zero(Tsubmatrix, dim)
            if is_det_zero:
                count_zero_det += 1
    
        if count_zero_det == len(submatrix):
            #get NEW submatrices (with lowest dimension) and do everything again
            new_submatrices = []
            for each_submatrix in submatrix:
                new_submatrices.extend(split_matrix(each_submatrix, dim-1)) # REVER!!
            return indirect_det(new_submatrices)
        else:
            return dim
    
    else:
        return 1


def solutionize(matrix, i, j):
    '''
    Indicates one of the three types of a solution in linear algebra system.
    
    matrix (list of list of floats): matrix organized by rows
    i (int): number of rows of the matrix 
    j (int): number of columns of the matrix 
    
    Returns: string with type of solution
    '''
    # first, find A
    A = []
    for k in range(i):
        A.append(matrix[k][:-1])
        
    # calculates the determinant of A and get the rank(A)
    submatrix_A = squares(A)
    p_A = indirect_det(submatrix_A)
    
    # calculates the determinant of M and get the rank(M)
    submatrix_M = squares(matrix)
    p_M = indirect_det(submatrix_M)
    
    # compare rank of A and M just once
    postos = (p_A == p_M)
    
    # identify solution type (REVER)
    if p_A == min(i,j-1): # j-1 = number of columns of matrix A
        if i == j-1:
            return "SPD"
        elif i < j-1:
            return "SPI"
        else: # i > j-1
            if postos:
                return "SPD"
            else:
                return "SI"
    else: #p_A < min{i,j-1}
        if postos:
            return "SPI"
        else:
            return "SI"


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
#M = 'matrix_SPD_2.txt'
#M = 'matrix_SPI.txt'
#M = 'matrix_SPI_2.txt'
#M = 'matrix_SI.txt'
#M = 'matrix_SI_2.txt'
#M = 'matrix_grande_cheia.txt'
#M = 'matrix_grande_vazia.txt'

test = Gauss(M)
test.get_initial_matrix()
print('\n', "Gauss result:", end="")
test.result()
#print(test.get_size())
