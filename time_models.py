#%%
import os
print(os.getcwd())
#%%
import numpy as np
from valve_models_pure_python.valve_model import sim_karnopp_pure_python
from valve_models_numba.numba_opt import sim_karnopp_numba
from valve_models_c_cython_wrap._valve_models_c_wrap_cython import valve_model_c_cython_wrap
import time

#%%
# dados.csv
# t,OP,P,x_kano,x_karnopp,x_lugre,x_gms
data = np.genfromtxt('./valve_models_c/test_data/dados.csv', delimiter=',')

# %%
Ts = 1e-3
dt = 1e-5
PARAM_VALVULA = np.array([1.6, 210490, 0.0445, 2550, 0, 0.029, 180000, 41368, 3.571428571428571], dtype=np.float64)
PARAM_ATRITO = np.array([700, 780, 125000, 5.0e-04], dtype=np.float64)
pos0 = np.array([0, 0], dtype=np.float64)

simdata = sim_karnopp_numba(data[0:1000,2], Ts, PARAM_VALVULA, PARAM_ATRITO,
    pos0, dt)

# %%

t1 = time.time()
simdata = valve_model_c_cython_wrap(data[:,2], Ts, PARAM_VALVULA, PARAM_ATRITO,
    pos0, dt)
t2 = time.time()
print(f'Simulation time: {t2-t1} s')

# %%
import matplotlib.pyplot as plt

plt.figure()
plt.plot(data[:,0], data[:,4], 'k')
plt.plot(simdata['t'], simdata['x'], 'r')
plt.legend(['matlab', 'python'])
plt.title('x')

plt.figure()
plt.plot(data[:,0], data[:,2], 'k')
plt.plot(simdata['t'], simdata['P'], 'r')
plt.legend(['matlab', 'python'])
plt.title('P')
plt.show()

#%%

print(f'Max x error {np.max(np.absolute(simdata["x"] - data[:,4]))}')
# %%
