# _________________________________________________________________
# SPDX-License-Identifier: MIT License
# For more information check at: https://spdx.org/licenses/MIT.html
# 
# Copyright (C) 2017
# Yuri Bastos Gabrich <yuribgabrich[at]gmail.com>
# _________________________________________________________________

def text2wolfram(input_text):
    '''
    Convert a matrix structured as text to WolframAlpha syntax.
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
