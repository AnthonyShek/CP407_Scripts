
# Browns compressor 4 stages

import numpy as np

gamma = 1.4
R = 54.29 # ft-lbf/()lbm-mol R)
g = 32.2 # ft/s2
mup = 0.48
Q=826.7786 # cfm
m = 60.92851491 # lb/min
T = (20.0+273.0)*9/5  # rankine
P = [1.0157* 14.504]  # psi
Ptarget = 30 * 14.504 # psi
eta_pre = 0.75 # polytropic efficiency
stages = 4
r = (Ptarget/P[0])**(1/stages) # pressure ratio

stg = 0
n= 1/(1-(gamma-1)/(gamma*eta_pre)) # polytropic exponent
Zavg = 1 # because ideal gas
Hp = Zavg*R*T * (1/(n-1)) * (r**(1-1/n)-1) # ft-lb/lb overall polytropic head
Hpp = Hp/stages
u2 = np.sqrt(Hpp*g/mup) # fps, impellar speed
d2 = 12 # inches, see fig 5-26, inlet vol 826 cfm
N = 60*12*u2/np.pi/d2
Qls = Q/((r**(1-1/stages))**(1/n))

flowcoeff_in = 700 * Q/N/(d2**3)
flowcoeff_ls = 700 * Qls/N/(d2**3)

eta_in = 0.71
eta_ls = 0.64
eta_avg = np.mean([eta_in,eta_ls])    
n= 1/(1-(gamma-1)/(gamma*eta_avg))
T2 = T*r**(1-1/n) *5/9 #K

Wp = m*Hp/(eta_avg*33000) # horsepower
Wp += 0.01*Wp
Wp = Wp * 0.7457 # kW