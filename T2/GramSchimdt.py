import copy as c

def inner_prod(v1, v2):
    '''
    Makes dot product between a pair of vectors.
    '''    
    result = 0
    for k in range(len(v1)):
        result += v1[k]*v2[k]
        
    return result


def proj(v2, v1):
    '''
    Makes the projection of v2 in v1.
    '''
    div = inner_prod(v2,v1) / inner_prod(v1,v1)
    v2_ = []
    for k in v2:
        v2_.append(div * k)
    return v2_

    
def orthog(v, ref):
    '''
    Makes the orthogonality between a vector and its reference (projection).
    '''
    if len(ref) == 0: # or v == ref == ref[0]
        return v
    # else
    v_ = c.deepcopy(v)
    for k in range(len(ref)-1, -1, -1):
        x = proj(v, orthog(ref[-1], ref[0:k]))
        
        for i in range(len(v)):
            v_[i] -= x[i]
    
    return v_
    

def GramSchimdt(beta):
    '''
    Makes Gram-Schimdt process.
    '''
    try:
        ort = [c.deepcopy(beta[0])]
        for k in range(len(beta)-1, 0, -1):
            ort.append(orthog(beta[k], beta[0:k]))
        
        # orthonormalising
        norm = c.deepcopy(ort)
        for j in norm:
            mod = (inner_prod(j,j))**0.5
            for i in range(len(j)):
                j[i] = round(j[i]/mod,2)
            
        return norm
    except:
        print("Something wrong with the size of your input vectors.")


#---------------------------------------
#columns
v1 = [2,5,8]
v2 = [3,6,9]
v3 = [4,7,10]

beta = [v1, v2, v3]
print(GramSchimdt(beta))
