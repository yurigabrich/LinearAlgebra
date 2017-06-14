def det(matrix):
    '''
    Calculate the determinant of a square matrix.
	
	Returns: the scalar related to that matrix.
    '''
    n = len(matrix)
    
    if n == 1:
    	return matrix

    if n == 2:
        return (matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0])
        
    #else --> recursive until n == 2 :: do cofactor method
    submatrices = combinations(matrix, n-1)
    
    result = 0
    for k in range(n):
        cofactor = (-1)**(k) * det(submatrices[k])
        result += (matrix[0][k] * cofactor)
    return result
