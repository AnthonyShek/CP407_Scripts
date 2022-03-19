def getMolFR(data):
    mol_FR = 0 # mol/s
    for molecule in data.keys():
        mol_FR += data[molecule]['n']
    return mol_FR

def getCp(data, T):
    T_sp = sp.symbols('T')
    mol_FR = getMolFR(data)
    Cp = 0
    for molecule in data_s.keys():
        data_s[molecule]['y'] = data_s[molecule]['n']/mol_FR
        
    for molecule in data.keys():
        data[molecule]['Cp_eqn'] = 0
        for n, coeff in enumerate(data[molecule]['Cp_coeff']):
            data[molecule]['Cp_eqn'] += coeff*T_sp**n
        tmp_np = sp.lambdify(T_sp, data[molecule]['Cp_eqn'])
        Cp += tmp_np(T)*data[molecule]['y']
    return Cp


# "Isothermal" Staged compressor

import numpy as np
import sympy as sp
import sympy.parsing.sympy_parser as spp

data_s = {
    'N2': {'n':15.77777778, 'Hf':0.0, 'MW':28.01}, # mol/s
    'O2': {'n':0.000128408, 'Hf':0.0, 'MW':32.00},
    'Ar': {'n':0.490429758,'Hf':0.0, 'MW':39.95},
    'H2O': {'n':5.39917E-05,'Hf':-242, 'MW':18.00}
    }

data_s['N2'] ['Cp_coeff'] = [31.150, -1.357e-2, 26.796e-6, -1.168e-8] # a, b, c, d cp in J/mol/K
data_s['O2']['Cp_coeff'] = [28.106, -3.68e-6, 17.459e-6, -1.065e-8]      
data_s['Ar'] ['Cp_coeff'] = [20.804, -3.211e-5, 51.665e-9, 0.0]
data_s['H2O'] ['Cp_coeff'] = [32.243, 19.238e-4, 10.555e-6, -3.596e-9]


gamma = 1.4
R = 8.3145
Q=[0.6]
T = [20.0+273.0]
P = [1.0157]
W=[]
power=[]
cooler_power=[]
Ptarget = 30
Tc = 133.63
Pc = 37.858
r = (Ptarget/P[0])**(1/4)

stg = 0
while P[stg]<Ptarget:
    
    #################
    Ep = spp.parse_expr(input(f'Stage:[{stg+1}] \nEnter Ep \nSee fig 3.6, at Q={Q[stg]:.3f} m3/s \n'))
    #Ep = 0.73 at Q = 3.627
    #Ep = 0.71 at Q = 1.858  
    #Ep = 0.68 at Q = 0.964  
    #Ep = 0.65 at Q = 0.510 
    #################
    
    m= (gamma-1)/gamma/Ep
    # if Ptarget/P[stg]>3.0:
    #     stg += 1
    #     P.append(3*P[stg-1])
    # elif Ptarget/P[stg] <= 3.0:
    #     stg += 1
    #     P.append(Ptarget)
    T.append( T[stg-1]*(P[stg]/P[stg-1])**m )
    Tr = (T[stg]+T[stg-1])/(2*Tc)
    Pr = (P[stg-1]+P[stg])/(2*Pc)
    
    #################
    Cp_corre = spp.parse_expr(input(f'Stage:[{stg}] \nEnter Cp correction \nSee fig 3.2, at Tr={Tr:.3f} Pr={Pr:.3f} \n'))
    #Cp_corre = 0.06 at Tr=2.782, Pr=0.054
    #Cp_corre = 0.13 at Tr=2.802, Pr=0.161 
    #Cp_corre = 0.39 at Tr=2.836, Pr=0.483 
    #Cp_corre = 0.90 at Tr=2.237, Pr=0.758 
    #################
    
    Cp_actual = getCp(data_s, (T[stg]+T[stg-1])/2) - Cp_corre
    
    ###############
    Z = spp.parse_expr(input(f'Stage:[{stg}] \nEnter Compressibility Z \nSee fig 3.8, at Tr={Tr:.3f} Pr={Pr:.3f} \n'))
    #Z = 1 at Tr=2.782, Pr=0.054
    #Z = 1 at Tr=2.802, Pr=0.161 
    #Z = 1 at Tr=2.836, Pr=0.483 
    #Z = 1 at Tr=2.237, Pr=0.758 
    ###############
    
    n = 1/(1-m)
    # polytrpoic compressor work
    W.append( Z*T[stg-1]*R * (n/(n-1)) * ((P[stg]/P[stg-1])**((n-1)/n)-1)) # J/mol
    # actual work -> power
    power.append( W[stg-1]/Ep *getMolFR(data_s)) # W
    # adjusting vol flowrate with ideal gas
    Q.append( Q[stg-1] /P[stg]*P[stg-1]  *  T[stg]/T[stg-1] )
    
    # intercooler duty
    T_sp = sp.symbols('T')
    deltaH=0
    for molecule in data_s.keys():
        tmp_np = sp.integrate(data_s[molecule]['Cp_eqn'], T_sp)
        tmp_np = sp.lambdify(T_sp, tmp_np)
        deltaH += data_s[molecule]['n']*(tmp_np(T[stg])-tmp_np(T[stg-1]))
        cooler_power.append(data_s[molecule]['n']*(tmp_np(T[stg])-tmp_np(T[stg-1])))
        
    print(f'\n     Stage:[{stg}] \n     Compressor Power:{power[stg-1]:.3f} W \n     Outlet Temperature:{T[stg]:.3f} K \n     Intercooler Duty:{deltaH:.3f} W')
    T[stg]=T[stg-1]
    
total_P=0
for i in range(len(power)):
    total_P += power[i]+cooler_power[i]
print(f'\nTotal system power consumption: {total_P:.3f} W')