
import numpy as np
import pylab as plt

t_dat = np.linspace(1,22,22)
n_dat = [148.214987736682,	295.365716246948,	441.461799239346,	586.512782025605,	730.528139141138,	873.517272797423,	1015.48951131542,	1156.45410819791,	1296.42024127855,	1435.39701233133,	1573.39344627339,	1710.41849065347,	1846.48101492309,	1981.58980943786,	2115.75358403302,	2248.98096596697,	2381.28049699829,	2512.66062930457,	2643.12971990679,	2772.69602317522,	2901.36768097453,	3029.15270996897]
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
plt.ylabel('Gas Moles Accumulated (mol)')