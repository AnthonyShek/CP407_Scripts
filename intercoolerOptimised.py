def get_cp(data):
    for molecule in data.keys():
        data[molecule]['Cp_eqn'] = 0
        for n, coeff in enumerate(data[molecule]['Cp_coeff']):
            data[molecule]['Cp_eqn'] += coeff*T**n
    return 0


import numpy as np
import sympy as sp

data_s = {
    'N2': {'n':56800, 'Hf':0.0, 'MW':28.01}, # mol/hr
    'O2': {'n':0.462267061, 'Hf':0.0, 'MW':32.00},
    'Ar': {'n':1765.547128,'Hf':0.0, 'MW':39.95},
    'H2O': {'n':0.194370139,'Hf':-242, 'MW':18.00}
    }

data_s['N2'] ['Cp_coeff'] = [31.150, -1.357e-2, 26.796e-6, -1.168e-8] # a, b, c, d cp in J/mol/K
data_s['O2']['Cp_coeff'] = [28.106, -3.68e-6, 17.459e-6, -1.065e-8]      
data_s['Ar'] ['Cp_coeff'] = [20.804, -3.211e-5, 51.665e-9, 0.0]
data_s['H2O'] ['Cp_coeff'] = [32.243, 19.238e-4, 10.555e-6, -3.596e-9]

T = sp.symbols('T')
Ts0, Ts1 = 425.055, 20+273 

deltaH = 0
get_cp(data_s)
for molecule in data_s.keys():
    tmp_np = sp.integrate(data_s[molecule]['Cp_eqn'], T)
    tmp_np = sp.lambdify(T, tmp_np)
    deltaH += data_s[molecule]['n']*(tmp_np(Ts1)-tmp_np(Ts0)) # J/hr
    
print(f'Cooler Duty: {deltaH/60**2/1000} kW')


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

Tlm = (Ts0-Tt1-Ts1+Tt0)/np.log((Ts0-Tt1)/(Ts1-Tt0))
R = (Ts0-Ts1)/(Tt1-Tt0)
S = (Tt1-Tt0)/(Ts0-Tt0)

Ft = np.sqrt(R*R+1)*np.log((1-S)/(1-R*S))/((R-1)*np.log((2-S*(R+1-np.sqrt(R*R+1)))/(2-S*(R+1+np.sqrt(R*R+1)))))
#Ft = 0.975 #see graph for two shell pass, four tube pass
Tm = Ft* Tlm

U = 60# W/m2/C see graph from coulson, find uncertainty
A = -deltaH/60**2/(Tm*U)

L = 1.83 #standard 12ft
Do = 16e-3 #standard steel tube
Di = Do-1.7e-3
At = L*np.pi*Do
Nt=A/At

K = 0.175 #using triangular pitch, 4 pass, higher tranfer rate
n = 2.285
D_bundle = Do*(Nt/K)**(1/n)
D_clear = 15e-3 #see fig 12.12 for shell bundle clearance, using fixed utube head, D_bundle = 0.658
D_shell = D_bundle + D_clear


##############################################################################
# TUBE SIDE COEFF
Tt_avg = np.mean([Tt0, Tt1])
CSA_t = np.pi/4*Di**2
Np = 4 # No. passes
tpp = Nt/Np # tubes per pass, 4 passes
flowA = tpp*CSA_t # total flow area

massv_t = data_t['H2O']['n']*18/1000/60**2/flowA # water mass velocity kg/s/m2
rhow = 997 # kg/m3
linv = massv_t/rhow # water linear velocity m/s

hi = 4200*(1.35+0.02*(Tt_avg-273))*linv**0.8/((Di*1000)**0.2)

##############################################################################
# SHELL SIDE COEFF
baff_space = D_shell*0.4*2 # 0.4 is optimal baffle spacing, factor of two to decrease pressure drop
t_pitch = 1.25*Do # pitch distance between tube centres
CSA_s = (t_pitch-Do)/Do*baff_space

massv_s = 0 # kg/s
for molecule in data_s.keys():
    massv_s += data_s[molecule]['n']*data_s[molecule]['MW']/(60**2)/1000
massv_s = massv_s/CSA_s

de = 1.1/Do*(t_pitch**2-0.917*Do**2)
Ts_avg = np.mean([Ts0, Ts1])

rho_air = 1.225 # kg/m3
mu_air = 2.18e-5 # Pa s avg betweem air temp range
Cp_air = 1 # kJ/kg/K
k_air = 31.85e-3 # W/m/K avg betweem air temp range

Re_s = massv_s * de/mu_air
Pr = Cp_air*1000 * mu_air/k_air
jh= 1.15e-2# using 25% baffle cut cos thats the optimum, see fig 12.29 Re = 2631
ho = k_air/de*jh*Re_s*Pr**(1/3)

Twall = Ts_avg - U/ho*(Ts_avg-Tt_avg)
mu_air_Tw = 19.07e-6 # visco at Twall=325.8K, see eng toolbox 

corr_factor = (mu_air/mu_air_Tw)**0.14
hs=ho*corr_factor

##############################################################################
# OVERALL COEFF
k_steel = 58 # W/m/K see coulson table 12.6
air_foul_coeff = 5000 # W/m2/C 
water_foul_coeff = 1000 # W/m2/C

U_actual= (1/hs + 1/air_foul_coeff + Do*np.log(Do/Di)/(2*k_steel) + (Do/Di)/water_foul_coeff + (Do/Di)/hi)**-1

##############################################################################
# P DROP TUBE
mu_w = 1e-3

Re_t = rhow*Di*linv/mu_w
jf_t = 6e-3 #see fig 12.24, Re = 4313
deltaPt = Np*(8*jf_t*(L/Di)+2.5)*rhow*linv**2/2

##############################################################################
# P DROP SHELL

linv_s = massv_s/rho_air
jf_s = 6e-2 #see fig 12.30, baffle cut 25% Re = 2631
deltaPs = 8*jf_s*(D_shell/de)*(L/baff_space)*rho_air*linv_s**2/2
