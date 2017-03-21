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
