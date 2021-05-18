#pragma once
#ifndef _VALVE_MODELS_H
#define _VALVE_MODELS_H


typedef struct results
{
    double* t;
    double* P;
    double* x;
    double* v;
    double* a;
    double* F_at;
    double* SP;
    double* OP;
    int len;
} results;

results* sim_karnopp(const double* P, const size_t len, double Ts, const double* param_valv, const double* param_atrito, const double* pos0, double dt);


#endif //_VALVE_MODELS_H