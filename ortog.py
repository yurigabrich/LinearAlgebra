def prod_int(v1, v2):
    '''
    Makes the usual product between a pair of vectors.
    '''
    #count the dimensions--------------------------------for now only linear, not matrix
    if type(v1) != list:
        v1 = [v1]
    if type(v2) != list:
        v2 = [v2]
    
    #test conditions
    if len(v1) != len(v2):
        return "Impossible, bro! Both vectors must have the same size."
    
    #make the product
    result = 0
    for k in range(len(v1)):
        result += v1[k]*v2[k]
        
    return result


def ortog(beta):
    '''
    Makes the ortogonomality of...
    '''
    if type(beta[-1]) == list:
        result = beta[-1][-1]
    else:
        result = beta[-1]
    
    for n in range(len(beta)-1, 0, -1):
        vn = beta[n]
        vn_1_ = ortog(beta[n-1])
        
        result = result - ( prod_int(vn, vn_1_) / prod_int(vn_1_, vn_1_) ) * vn_1_
        
    return result

%-----------------------------------------------    
v1 = [2]
v2 = [3]
v3 = [4]

beta = [v1, v2, v3]



https://goo.gl/wQuTCU
