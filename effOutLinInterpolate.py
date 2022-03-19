
import numpy as np
import pylab as plt

t_dat = np.linspace(1,22,22)
n_dat = [18.8050122633185, 19.8692714897333,20.9239170076018,21.9690172137412,23.0046428844673,24.0308663437148,25.0477614820005,26.0554031175138,27.0538669193602,28.0432289472174,29.0235660579368,29.9949556199268,30.9574757303768,31.9112054852264,32.8562254048477,33.7926180660474,34.7204689686776,35.6398676937250,36.5509093977815,37.4536967315667,38.3483422006903,39.2349710055598]
plt.plot(t_dat, n_dat, 'x')

t_avg, n_avg = np.mean(t_dat), np.mean(n_dat)

nume, denom = 0,0
for t, n in zip(t_dat, n_dat):
    nume += (t-t_avg)*(n-n_avg)
    denom += (t-t_avg)**2
slope= nume/denom
intercept = n_avg - slope*t_avg
print(f'y = {slope:.4f}x+{intercept:.4f}')

t_plt = np.linspace(min(t_dat), max(t_dat), 10)
n_plt = [(intercept+slope*t) for t in t_plt]

plt.plot(t_plt, n_plt, label='Least Square Line')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Effluent Molar Flowrate (mol/s)')