def list2wolfram(input_list):
    '''
    Convert a matrix structured as list of list to WolframAlpha syntax.
    '''
    
    output = '{{'
    for i in input_list:
        for j in i:
            output += str(j)
            if i.index(j) != len(i)-1:
                output += ","
        if i != input_list[-1]:
            output += "},{"
    output += '}}'
    
    return output
    
input_list = [[1.0, 6.0, 10.0], [2.0, 7.0, 11.0], [3.0, 8.0, 12.0]]

print(list2wolfram(input_list))
