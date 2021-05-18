
import numpy as np
from typing import Dict


def initalize_simdata(simdata: Dict[str, np.ndarray], len: int) -> np.ndarray:
    simdata['t'] = np.zeros(len, dtype=np.float64)
    simdata['x'] = np.zeros(len, dtype=np.float64)
    simdata['P'] = np.zeros(len, dtype=np.float64)
    simdata['v'] = np.zeros(len, dtype=np.float64)
    simdata['a'] = np.zeros(len, dtype=np.float64)
    simdata['F_at'] = np.zeros(len, dtype=np.float64)
    simdata['SP'] = np.zeros(len, dtype=np.float64)
    simdata['OP'] = np.zeros(len, dtype=np.float64)
    return simdata


def upsample(x_init: float, x_end: float, nSamp: int, return_data: np.ndarray) -> np.ndarray:
    for i in range(nSamp):
        return_data[i] = x_init + (x_end - x_init) * i / nSamp

    return return_data

def sim_karnopp_pure_python(u: np.ndarray, Ts: float, 
                            param_valv: np.ndarray, param_atrito: np.ndarray, 
                            pos0: np.ndarray, dt: float) -> Dict[str, np.ndarray]:


    m = param_valv[0]
    k = param_valv[1]
    S_a = param_valv[2]
    F_init = param_valv[3]
    x_min = param_valv[4]
    x_max = param_valv[5]
    p_max = param_valv[6]
    p_min = param_valv[7]
    tau_ip = param_valv[8]

    F_c = param_atrito[0]
    F_s = param_atrito[1]
    F_v = param_atrito[2]
    v_s = param_atrito[3]

    nSamp = int(Ts/dt)
    len = u.shape[0]

    x = np.zeros(nSamp+2, dtype=np.float64)
    x[nSamp+1] = pos0[0]
    x[nSamp] = pos0[0]

    OP_exc = np.zeros(nSamp+2, dtype=np.float64)
    P = np.zeros(nSamp+2, dtype=np.float64)
    P_exc = np.zeros(nSamp+2, dtype=np.float64)
    P_us = np.zeros(nSamp, dtype=np.float64)
    v = np.zeros(nSamp+2, dtype=np.float64)
    a = np.zeros(nSamp+2, dtype=np.float64)
    F_at = np.zeros(nSamp+2, dtype=np.float64)
    F_res = np.zeros(nSamp+2, dtype=np.float64)
    
    simdata = dict()
    simdata = initalize_simdata(simdata, len)

    stick = 1
    P_exc_old = 0
    P_exc_old_old = 0

    for i in range(len):

        if i == len - 1:
            simdata['x'][i] = x[nSamp+1]
            simdata['v'][i] = v[nSamp+1]
            simdata['a'][i] = a[nSamp+1]
            simdata['P'][i] = u[len-1]
            simdata['F_at'][i] = F_at[nSamp+1]
            simdata['t'][i] = i*Ts
            break


        P_us = upsample(u[i], u[i+1], nSamp, P_us)
        P_exc_old = P_exc[nSamp+1]
        P_exc_old_old = P_exc[nSamp]
        if (i == 0):
            P_exc[0] = P[0]
            P_exc[1] = P[0]
        else:
            P_exc[0] = P_exc_old_old
            P_exc[1] = P_exc_old
		
        for ct in range(nSamp):
            P_exc[ct+2] = P_us[ct]
		
        x[0] = x[nSamp]
        x[1] = x[nSamp+1]
        v[0] = v[nSamp]
        v[1] = v[nSamp+1]
        a[0] = a[nSamp]
        a[1] = a[nSamp+1]
        F_at[0] = F_at[nSamp]
        F_at[1] = F_at[nSamp+1]
        F_res[0] = F_res[nSamp]
        F_res[1] = F_res[nSamp+1]

        simdata['x'][i] = x[1]
        simdata['v'][i] = v[1]
        simdata['a'][i] = a[1]
        simdata['P'][i] = P_exc[2]
        simdata['F_at'][i] = F_at[1]
        simdata['t'][i] = float(i) * Ts

        if np.isnan(simdata['x'][i]):
            simdata = initalize_simdata(simdata, len)
            return simdata

        for j in range(nSamp+2):
            if v[j-2] == 0 and not(stick):
                stick = 0
            elif v[j-2]*v[j-1] <= 0:
                stick = 1
            else:
                stick = 0

            if stick:
                F_r = S_a*P_exc[j-1] - k*x[j-1] - F_init
                if F_r > 0:
                    sig_F = 1
                elif F_r < 0:
                    sig_F = -1
                else:
                    sig_F = 0
                F_at[j] = sig_F*np.minimum(np.absolute(F_r), F_s)
                v[j-1] = 0
                if np.absolute(F_at[j]) >= F_s:
                    stick = 0
            else:
                if v[j-1] > 0:
                    sinal = 1
                elif v[j-1] < 0:
                    sinal = -1
                else:
                    sinal = 0
                F_at[j] = ( F_c + (F_s - F_c)*np.exp( -np.power(v[j-1]/v_s, 2) ) )*sinal + F_v*v[j-1]

            F_res[j] = S_a*P_exc[j-1] - k*x[j-1] - F_init - F_at[j]

            a[j] = 1/m*F_res[j]
            v[j] = dt/2 * (a[j] + a[j-1]) + v[j-1]
            x[j] = dt/2 * (v[j] + v[j-1]) + x[j-1]

            if x[j] < x_min:
                F_res[j] = 0
                a[j] = 0
                v[j] = 0
                x[j] = x_min
                stick = 1
            elif x[j] > x_max:
                F_res[j] = 0
                a[j] = 0
                v[j] = 0
                x[j] = x_max
                stick = 1

    return simdata

