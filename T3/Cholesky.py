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
        for n in range(i_turn, len(aug_matrix)): # because 'i' changes recursively!
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
                # change the value of entirely row
                num = aug_matrix[i_turn][k]
                denom = aug_matrix[k][k]
                for x in range(len(aug_matrix[0])): # because 'j' changes recursively!
                    aug_matrix[i_turn][x] -= (aug_matrix[k][x] * (num / denom))
            # UPDATE PIVOT
            pivot = aug_matrix[i_turn][j_turn]
        
        #stop criterion
        if min(i,j) < 2:
            stop = True
        
        return echelon(aug_matrix, i-1, j-1, stop)
    
    else:
        return aug_matrix #stepped
        

def is_square(matrix):
    '''
    Determines if a matrix has the same number of rows and columns.
    
    Returns: True if is a square matrix, False otherwise.
    '''
    return (len(matrix) == len(matrix[0]))
    

def is_symmetric(matrix):
    '''
    Determines if a matrix is symmetric, that is (matrix) = (matrix)^T
    
    Returns: True if is a symmetric matrix, False otherwise.
    '''
    if is_square:
        n = len(matrix)
        for i in range(n):
            for j in range(i, n):
                if matrix[i][j] != matrix[j][i]:
                    return False
        return True
    # else:
    return False
    

def is_positive(matrix):
    '''
    Determines if a matrix is definite-positive analysing its pivots.
    
    Returns: True if is a definite-positive matrix, False otherwise.
    '''
    reduced = echelon(matrix, len(matrix[0]), len(matrix))
    
    for row in reduced:
        for pivot in row:
            if pivot != 0:
                # it's definitely a pivot
                if pivot < 0:
                    return False
                break
                # Will it stop only 1 loop and go straight to the next row?!
    return True


def Cholesky(matrix):
    '''
    Calculates the Cholesky factor (G).
	
	Returns: the matrix G (list of lists).
    '''
    try:
        if is_positive(matrix) and is_symmetric(matrix):
            A = matrix
            n = len(A)
            G = []
            
            for k in range(n):
                # Firstly, create G
                G.append([])
                for j in range(n):
                    G[k].append(0)
                
                # Then, populate it
                sum = 0
                for i in range(k-1):
                    sum += G[k][i]**2
                    
                G[k][k] = round( (A[k][k] - sum)**0.5 ,2)
                
            for k in range(n):
                for j in range(k+1, n):
                    sum = 0
                    for i in range(k-1):
                        sum += G[j][i] * G[k][i]
                        
                    G[j][k] = round( (1/G[k][k]) * (A[j][k] - sum) ,2)
            
            return G
    
    except:
        return "The Cholesky factor couldn't be computed. Check inputted matrix."


#---------------------------------------
#rows
v1 = [2,5,8]
v2 = [3,6,9]
v3 = [4,7,10]

matrix = [v1, v2, v3]
x = Cholesky(matrix)
if type(x)!=str():
  for i in x: print(i)
else: print(x)
