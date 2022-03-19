def get_cp(data):
    for molecule in data.keys():
        data[molecule]['Cp_eqn'] = 0
        for n, coeff in enumerate(data[molecule]['Cp_coeff']):
            data[molecule]['Cp_eqn'] += coeff*T**n
    return 0

import sympy as sp


deltaH = -3758 *1000*60**2 # J/hr mean duty
T = sp.symbols('T')
data_t = {
    'H2O': {'n':0.0,'Hf':-242, 'MW':18.00}
    }
data_t['H2O'] ['Cp_coeff'] = [32.243, 19.238e-4, 10.555e-6, -3.596e-9]
get_cp(data_t)
Tt0, Tt1 = 13.5+273, 13.5+8+273
for molecule in data_t.keys():
    tmp_np = sp.integrate(data_t[molecule]['Cp_eqn'], T)
    tmp_np = sp.lambdify(T, tmp_np)
data_t['H2O']['n'] = -deltaH/(tmp_np(Tt1)-tmp_np(Tt0)) # mol/hr 
data_t['H2O']['m'] = data_t['H2O']['n'] * data_t['H2O']['MW'] /1000/1000 *24 *365 # ton/yr