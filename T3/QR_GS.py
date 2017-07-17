def QR_GS(A):
    '''
    Calculates the QR decomposition by Gram-Schimdt method.
	
	Returns: Q, R (tuple of lists of list of floats)
    '''
    # calculating Q
    Q = GramSchimdt(A)
    
    # calculating R
    R = []
    for j in range(len(A)):
        R.append([])
        for i in range(len(A[0])):
            if i <= j:
               value = round(inner_prod(A[j], Q[i]), 2)
            else:
                value = 0
            R[j].append(value)
            
    return Q, R
