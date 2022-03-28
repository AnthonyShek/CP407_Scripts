
def getDiff(func, var):
    return sp.diff(func, var)

import sympy as sp
import numpy as np

p, s = sp.symbols('p s') #production rate, sale price
variables = [p, s]
p_dat, s_dat = 18369, 990
errors = [0.2*p_dat, 150] 
r = p*s #revenue

diffs = []
for var in variables:
    diffs.append(getDiff(r, var))

err_prop=0
for diff, error in zip(diffs, errors):
    err_prop += diff**2*error**2
    
r = sp.lambdify(variables, r)    
err_prop = sp.lambdify(variables, sp.sqrt(err_prop))

print(f'revenue: {r(p_dat, s_dat):.3f}±{err_prop(p_dat, s_dat):.3f}')