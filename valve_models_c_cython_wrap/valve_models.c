#include "stdio.h"
#include "stdlib.h"
#include "errno.h"
#include "math.h"
#include "valve_models.h"
#include "time.h"
#include "std_random.h"

void upsample(double x_init, double x_end, size_t nSamp, double *return_data);
void clear_simdata(results *simdata, size_t len);
void allocate_results(results *simdata, size_t len);
double nrand();
void 
calc_F_at(const double *friction_params, 
          const double *valve_params, 
          const char* friction_model, 
          const double* P_exc, 
          const double* x, 
          double* v, 
          double* F_at, 
          int* stick, 
          const size_t j);

void allocate_results(results *simdata, size_t len) {
	simdata->t = calloc(len, sizeof(double));
	simdata->x = calloc(len, sizeof(double));
	simdata->P = calloc(len, sizeof(double));
	simdata->v = calloc(len, sizeof(double));
	simdata->a = calloc(len, sizeof(double));
	simdata->F_at = calloc(len, sizeof(double));
    simdata->SP = calloc(len, sizeof(double));
    simdata->OP = calloc(len, sizeof(double));
	simdata->len = len;
}

void clear_simdata(results *simdata, size_t len) {
    for (size_t i = 0; i < len; i++) {
        simdata->t[i] = 0.0;
        simdata->x[i] = 0.0;
        simdata->P[i] = 0.0;
        simdata->v[i] = 0.0;
        simdata->a[i] = 0.0;
        simdata->F_at[i] = 0.0;
        simdata->SP[i] = 0.0;
        simdata->OP[i] = 0.0;
    }
}

double nrand() {
    double urand = ( (double) random() ) / ( (double) RAND_MAX );
    return ltqnorm(urand);
}


results* sim_karnopp(const double* u, const size_t len, double Ts, const double* param_valv, const double* param_atrito, const double* pos0, double dt) {
    // dt = simulation integration time
	// Ts = data sampling time
    if (Ts < dt){
		perror("Error: simulation sampling time has to be lower than data sampling time");
	}
        
	double m = param_valv[0];
	double k = param_valv[1];
	double S_a = param_valv[2];
	double F_init = param_valv[3];
	double x_min = param_valv[4];
	double x_max = param_valv[5];
	double p_max = param_valv[6];
	double p_min = param_valv[7];
	double tau_ip = param_valv[8];
    
	double F_c = param_atrito[0];
	double F_s = param_atrito[1];
	double F_v = param_atrito[2];
	double v_s = param_atrito[3];

	size_t nSamp = Ts/dt;

	// Valve model with Karnopp friction model
    // valve stem position initialization
	double *x = calloc(nSamp+2, sizeof(double));
	x[nSamp+1] = pos0[0];
	x[nSamp] = pos0[0];

    // OP initialization
    double* OP_exc = calloc(nSamp+2, sizeof(double));
	double* P = calloc(nSamp+2, sizeof(double));

	// diaphragm pressure initialization
	double *P_exc = calloc(nSamp+2, sizeof(double));
	double *P_us = calloc(nSamp, sizeof(double));
    
	// valve stem velocity initialization
	double *v = calloc(nSamp+2, sizeof(double));
	
	// valve stem acceleration initialization
	double *a = calloc(nSamp+2, sizeof(double));
	
    // friction force initialization
	double *F_at = calloc(nSamp+2, sizeof(double));

	// resultant force initialization
	double *F_res = calloc(nSamp+2, sizeof(double));

	results *simdata = malloc(sizeof(results));
    allocate_results(simdata, len);
	
    int stick = 1;
	double P_exc_old, P_exc_old_old;

	for (size_t i = 0; i < len; i++){

		if (i == len-1){
			simdata->x[i] = x[nSamp+1];
			simdata->v[i] = v[nSamp+1];
			simdata->a[i] = a[nSamp+1];
			simdata->P[i] = u[len-1];
			simdata->F_at[i] = F_at[nSamp+1];
			simdata->t[i] = i*Ts;
			break;
		}

		// reinitialize intermediate vectors
		upsample(u[i], u[i+1], nSamp, &P_us[0]);
		P_exc_old = P_exc[nSamp+1];
		P_exc_old_old = P_exc[nSamp];
		if (i == 0){
			P_exc[0] = P[0];
			P_exc[1] = P[0];
		}
		else{
			P_exc[0] = P_exc_old_old;
			P_exc[1] = P_exc_old;
		}
		for (int ct = 0; ct < nSamp; ct++){
			P_exc[ct+2] = P_us[ct];
		}

		
		x[0] = x[nSamp];
		x[1] = x[nSamp+1];
		v[0] = v[nSamp];
		v[1] = v[nSamp+1];
		a[0] = a[nSamp];
		a[1] = a[nSamp+1];
		F_at[0] = F_at[nSamp];
		F_at[1] = F_at[nSamp+1];
		F_res[0] = F_res[nSamp];
		F_res[1] = F_res[nSamp+1];

		simdata->x[i] = x[1];
		simdata->v[i] = v[1];
		simdata->a[i] = a[1];
		simdata->P[i] = P_exc[2];
		simdata->F_at[i] = F_at[1];
		simdata->t[i] = (double) (i*Ts);

        if (isnan(simdata->x[i])) {
			free(x);
			free(P_exc);
			free(P_us);
			free(v);
			free(a);
			free(F_at);
			free(F_res);
			free(OP_exc);
			free(P);
			clear_simdata(simdata, len);
			return simdata;
		}

		for (size_t j = 2; j < nSamp+2; j++){
 
           calc_F_at(param_valv, param_atrito, "teste", &P_exc[0], &x[0], &v[0], &F_at[0], &stick, j);

            F_res[j] = S_a*P_exc[j-1] - k*x[j-1] - F_init - F_at[j];

			a[j] = 1/m*F_res[j];
			v[j] = dt/2 * (a[j] + a[j-1]) + v[j-1];
			x[j] = dt/2 * (v[j] + v[j-1]) + x[j-1];

			if (x[j] < x_min){
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_min;
                stick = 1;
			}
			else if (x[j] > x_max){
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_max;
                stick = 1;
			}
		}

	}

    free(x);
    free(P_exc);
    free(P_us);
    free(v);
    free(a);
    free(F_at);
    free(F_res);
    free(OP_exc);
    free(P);
	return simdata;
}

void 
calc_F_at(const double *param_valv, 
          const double *param_atrito,
          const char* friction_model, 
          const double* P_exc, 
          const double* x, 
          double* v, 
          double* F_at, 
          int* stick,
          const size_t j) {

	double m = param_valv[0];
	double k = param_valv[1];
	double S_a = param_valv[2];
	double F_init = param_valv[3];
	double x_min = param_valv[4];
	double x_max = param_valv[5];
	double p_max = param_valv[6];
	double p_min = param_valv[7];
	double tau_ip = param_valv[8];
    
	double F_c = param_atrito[0];
	double F_s = param_atrito[1];
	double F_v = param_atrito[2];
	double v_s = param_atrito[3];

    double sinal, F_r, sig_F;
    if (v[j-2] == 0 && !*stick)
        *stick = 0;
    else if (v[j-2]*v[j-1] <= 0)
        *stick = 1;
    else
        *stick = 0;

    if (*stick){
        F_r = S_a*P_exc[j-1] - k*x[j-1] - F_init;
        if (F_r > 0)
            sig_F = 1;
        else if (F_r < 0)
            sig_F = -1;
        else
            sig_F = 0;
        F_at[j] = sig_F*fmin(fabs(F_r), F_s);
        v[j-1] = 0;
        if (fabs(F_at[j]) >= F_s)
            *stick = 0;
    }
    else{
        if (v[j-1] > 0)
            sinal = 1;
        else if (v[j-1] < 0)
            sinal = -1;
        else
            sinal = 0;
        F_at[j] = ( F_c + (F_s - F_c)*exp( -pow(v[j-1]/v_s, 2) ) )*sinal + F_v*v[j-1];
    }
}

void upsample(double x_init, double x_end, size_t nSamp, double *return_data){
	for (int i = 0; i < nSamp; i++){
		return_data[i] = x_init + (x_end - x_init) * i / nSamp;
	}
}
