def text2wolfram(input_text):
    '''
    Convert a matrix structure as text to
    WolframAlpha sintax.
    '''
    
    output = '{{'
    for i in input_text:
        if i == ' ':
            output += ','
        elif i == ',':
            output += '},{'
        else:
            output += i
    output += '}}'
    
    return output
    
input_text = "2 1 3,1 -3 -2"

print(text2wolfram(input_text))
