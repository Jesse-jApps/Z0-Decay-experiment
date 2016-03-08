import numpy as np 
import re 

def _latex_number(number, prec=10):
    if number == 0:
        return '0'

    stringified = ('{:.'+str(prec)+'e}').format(number)
    if 'e' in stringified:
        match = re.match(u'([+-]?\d*\.?\d*)e([+-]\d+)', stringified)
        
        base = ''
        if float(match.group(2)) != 0:
            base = '\cdot 10^{'+match.group(2)+'}'
        
        return str(round(float(match.group(1)),prec))+base

    return str(round(number,prec))

def np_array_to_latex(data, precision=1000):
    ltx = ['\\begin{pmatrix}']
    ltx += [('   ' + ' & '.join(_latex_number(a,precision) for a in r)+' \\\\') for r in data]
    ltx += ['\\end{pmatrix}']
    return '\n'.join(ltx)

if __name__ == '__main__':
    data = np.array([[1,2,3],[3.3,4,5], [4,5,6]])
    print(np_array_to_latex(data))