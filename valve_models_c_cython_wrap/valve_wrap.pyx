
import numpy as np
cimport numpy as np
np.import_array()


cdef extern from "valve_models.h":
    ctypedef struct results:
        double* t
        double* P
        double* x
        double* v
        double* a
        double* F_at
        double* SP
        double* OP
        int len
    results* sim_karnopp(const double* P, const size_t len, double Ts, const double* param_valv, const double* param_atrito, const double* pos0, double dt)


def valve_model_c_cython_wrap(np.ndarray[double, ndim=1] u, 
                              double Ts, 
                              np.ndarray[double, ndim=1] param_valv,
                              np.ndarray[double, ndim=1] param_atrito, 
                              np.ndarray[double, ndim=1] pos0, 
                              double dt):
    
    u = np.ascontiguousarray(u)
    param_valv = np.ascontiguousarray(param_valv)
    param_atrito = np.ascontiguousarray(param_atrito)
    pos0 = np.ascontiguousarray(pos0)

    cdef double[::1] u_memview = u
    cdef double[::1] param_valv_memview = param_valv
    cdef double[::1] param_atrito_memview = param_atrito
    cdef double[::1] pos0_memview = pos0
    cdef results* simdata
    simdata = sim_karnopp(&u_memview[0], 
                          u_memview.shape[0], 
                          Ts, 
                          &param_valv_memview[0], 
                          &param_atrito[0], 
                          &pos0[0], 
                          dt)


    retval = {
        't': np.zeros(u_memview.shape[0], dtype=np.float64),
        'P': np.zeros(u_memview.shape[0], dtype=np.float64),
        'x': np.zeros(u_memview.shape[0], dtype=np.float64),
        # 'v': np.zeros(u_memview.shape[0], dtype=np.float64),
        # 'a': np.zeros(u_memview.shape[0], dtype=np.float64),
        # 'F_at': np.zeros(u_memview.shape[0], dtype=np.float64),
        'SP': np.zeros(u_memview.shape[0], dtype=np.float64),
        'OP': np.zeros(u_memview.shape[0], dtype=np.float64),
    }
    cdef int i
    for i in range(u_memview.shape[0]):
        retval['t'][i] = simdata.t[i]
        retval['P'][i] = simdata.P[i]
        retval['x'][i] = simdata.x[i]
        # retval['v'][i] = simdata.v[i]
        # retval['a'][i] = simdata.a[i]
        # retval['F_at'][i] = simdata.F_at[i]
        retval['SP'][i] = simdata.SP[i]
        retval['OP'][i] = simdata.OP[i]
    return retval