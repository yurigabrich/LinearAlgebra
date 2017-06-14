import copy as c

def inner_prod(v1, v2):
    '''
    Makes dot product between a pair of vectors.
    '''    
    result = 0
    for k in range(len(v1)):
        result += v1[k]*v2[k]
        
    return result


def ortog(beta):
    '''
    Makes the orthogonality base of a set of vectors.
    '''
    if type(beta[0]) != list:
        result = c.deepcopy(beta) #v1_ = v1
    else:
        result = c.deepcopy(beta[-1]) #vn_
    
        for n in range(len(beta)-1, 0, -1):
            vn = beta[n]
            vn_1_ = ortog(beta[n-1])
            frac = inner_prod(vn, vn_1_) / inner_prod(vn_1_, vn_1_)
            
            for k in range(len(result)):
                result[k] -= (frac * vn_1_[k])
                result[k] = round(result[k], 2)
        
    return result
    

def GramSchimdt(beta):
    '''
    Makes Gram-Schimdt process.
    '''
    Try:
        return ortog(beta) #still missing  orthonormalising a set of vectors
    Except:
        print("Something wrong with the size of your input vectors.")


%---------------------------------------
v1 = [2,5,8]
v2 = [3,6,9]
v3 = [4,7,10]

beta = [v1, v2, v3]
print(ortog(beta))

https://goo.gl/ONDBXL
