results* sim_karnopp(const double* u, const size_t len, double Ts, const double* param_valv, const double* param_atrito, const double* pos0, double dt, double t0) {
	// dt = simulation integration time
	// Ts = data sampling time

    if (Ts < dt)
		perror("Error: simulation sampling time has to be lower than data sampling time");

    char* simulation_type = "ol";
            
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

	size_t nSamp = (size_t) (Ts/dt);

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
	results* simdata = malloc(sizeof(results));

    allocate_results(simdata, len);

    bool stick = true;
	double F_r = 0, sig_F = 0, P_exc_old = 0, P_exc_old_old = 0;
    double Fs_m_Fc = F_s - F_c, one_m = 1 / m, dt_2 = dt / 2;

    srandom(time(NULL));

	for (size_t i = 0; i < len; i++){

        if (i == len - 1) {
			simdata->x[i] = x[nSamp + 1];
            simdata->v[i] = v[nSamp + 1];
            simdata->a[i] = a[nSamp + 1];
            simdata->F_at[i] = F_at[nSamp +1];
            simdata->t[i] = i * Ts + t0;
			if (!strcmp(simulation_type, "ol")) {
				simdata->P[i] = u[len - 1];
			}
			else if (!strcmp(simulation_type, "cl")) {
				simdata->P[i] = P_exc[nSamp + 1];
				simdata->OP[i] = OP_exc[nSamp + 1];
				simdata->SP[i] = u[i];
			}
			break;
		}

		// reinitialize intermediate vectors
		if (!strcmp(simulation_type, "ol")) {
			upsample(u[i], u[i + 1], nSamp, &P_us[0]);
			P_exc_old = P_exc[nSamp + 1];
			P_exc_old_old = P_exc[nSamp];
            if (i == 0) {
				P_exc[0] = u[0];
				P_exc[1] = u[0];
			}
			else {
				P_exc[0] = P_exc_old_old;
				P_exc[1] = P_exc_old;
			}
			for (size_t ct = 0; ct < nSamp; ct++) {
				P_exc[ct + 2] = u[ct];
			}
		}
		else if (!strcmp(simulation_type, "cl")) {
			double OP = 0;
			if (i == 0) {
				OP_exc[0] = 0;
				OP_exc[1] = 0;
			}
			else {
				OP = 100;//controller.pid(u[i], (simulation_results.x[i - 1] - x_min) / (x_max - x_min) * 100 + nrand(), i);
				if (OP > 100)
					OP = 100;
				else if (OP < 0)
					OP = 0;
				OP_exc[0] = simdata->OP[i - 1];
				OP_exc[1] = simdata->OP[i - 1];
			}
			for (size_t ct = 2; ct < nSamp + 2; ++ct) {
				OP_exc[ct] = OP;
				P[ct] = (dt * (p_max - p_min) / 100 * OP_exc[ct - 1] + tau_ip * P[ct - 1]) / (dt + tau_ip);
				P_exc[ct] = P[ct] + p_min;
			}
			P[0] = P[nSamp];
			P[1] = P[nSamp + 1];
			P_exc[0] = P_exc[nSamp];
			P_exc[1] = P_exc[nSamp + 1];
			simdata->OP[i] = OP;
			simdata->SP[i] = u[i];
		}

		


 

		for (size_t j = 2; j < nSamp + 2; j++) {
			if (v[j - 2] == 0 && !stick)
				stick = false;
			else if (v[j - 2] * v[j - 1] <= 0)
				stick = true;
			else
				stick = false;
			if (stick) {
				F_r = S_a * P_exc[j - 1] - k * x[j - 1] - F_init;
				double absFr = (F_r < 0 ? -F_r : F_r);
				if (F_r > 0)
					F_at[j] = (absFr < F_s ? absFr : F_s);
				else if (F_r < 0)
					F_at[j] = (absFr < F_s ? -absFr : -F_s);
				else
					F_at[j] = 0;

				F_res[j] = F_r - F_at[j];
				v[j - 1] = 0;

				double absFat = (F_at[j] < 0 ? -F_at[j] : F_at[j]);
				stick = (absFat >= F_s ? false : true);
			}
			else {
				if (v[j - 1] > 0) {
					double v_vs = v[j - 1] / v_s;
					double exp_v2 = exp(-v_vs * v_vs);
					F_at[j] = (F_c + Fs_m_Fc * exp_v2) + F_v * v[j - 1];
				}
				else if (v[j - 1] < 0) {
					double v_vs = v[j - 1] / v_s;
					double exp_v2 = exp(-v_vs * v_vs);
					F_at[j] = -(F_c + Fs_m_Fc * exp_v2) + F_v * v[j - 1];
				}
				else {
					F_at[j] = 0;
				}
				F_res[j] = S_a * P_exc[j - 1] - k * x[j - 1] - F_at[j] - F_init;
			}

			a[j] = one_m * F_res[j];
			v[j] = dt_2 * (a[j] + a[j - 1]) + v[j - 1];
			x[j] = dt_2 * (v[j] + v[j - 1]) + x[j - 1];

			if (x[j] < x_min) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_min;
			}
			else if (x[j] > x_max) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_max;
			}
		}
	}

    // Free the variables
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





            printf("%e, %e\n", u[i], u[i+1]);
            printf("%e, %e\n", P_us[0], P_us[nSamp-1]);
            printf("%d\n", nSamp);
            printf("----------------------------\n");



	for (size_t i = 0; i < u.size(); i++) {


	


		if (isnan(simulation_results.x[i])) {
			delete[] x;
			delete[] P_exc;
			delete[] P_us;
			delete[] v;
			delete[] a;
			delete[] F_at;
			delete[] F_res;
			delete[] OP_exc;
			delete[] P;
			clear_sim_data();
			return;
		}

		for (size_t j = 2; j < nSamp + 2; j++) {
			if (v[j - 2] == 0 && !stick)
				stick = false;
			else if (v[j - 2] * v[j - 1] <= 0)
				stick = true;
			else
				stick = false;
			if (stick) {
				F_r = S_a * P_exc[j - 1] - k * x[j - 1] - F_init;
				real absFr = (F_r < 0 ? -F_r : F_r);
				if (F_r > 0)
					F_at[j] = (absFr < F_s ? absFr : F_s);
				else if (F_r < 0)
					F_at[j] = (absFr < F_s ? -absFr : -F_s);
				else
					F_at[j] = 0;

				F_res[j] = F_r - F_at[j];
				v[j - 1] = 0;

				real absFat = (F_at[j] < 0 ? -F_at[j] : F_at[j]);
				stick = (absFat >= F_s ? false : true);
			}
			else {
				if (v[j - 1] > 0) {
					real v_vs = v[j - 1] / v_s;
					real exp_v2 = std::exp(-v_vs * v_vs);
					F_at[j] = (F_c + Fs_m_Fc * exp_v2) + F_v * v[j - 1];
				}
				else if (v[j - 1] < 0) {
					real v_vs = v[j - 1] / v_s;
					real exp_v2 = std::exp(-v_vs * v_vs);
					F_at[j] = -(F_c + Fs_m_Fc * exp_v2) + F_v * v[j - 1];
				}
				else {
					F_at[j] = 0;
				}
				F_res[j] = S_a * P_exc[j - 1] - k * x[j - 1] - F_at[j] - F_init;
			}

			a[j] = one_m * F_res[j];
			v[j] = dt_2 * (a[j] + a[j - 1]) + v[j - 1];
			x[j] = dt_2 * (v[j] + v[j - 1]) + x[j - 1];

			if (x[j] < x_min) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_min;
			}
			else if (x[j] > x_max) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_max;
			}
		}

	}
	delete[] x;
	delete[] P_exc;
	delete[] P_us;
	delete[] v;
	delete[] a;
	delete[] F_at;
	delete[] F_res;
	delete[] OP_exc;
	delete[] P;
}



template <typename real>
void ValveModel<real>::sim_kano() {

	real* du = new real[u.size()]{ 0 };
	real* x = new real[u.size()]{ 0 };
	real* input = new real[u.size()]{ 0 };
	real* input_int = new real[u.size()]{ 0 };
	real* noise = new real[u.size()]{ 0 };
	real* OP = new real[u.size()]{ 0 };
	real x0 = pos0[0];

	allocate_sim_data(u.size());

	int stp = 1;
	input[0] = (u[0] - p_min) / (p_max - p_min) * 100;
	input[0] = (input[0] < 0) ? 0 : input[0];
	input[0] = (input[0] > 100) ? 100 : input[0];
	input[0] = input[0] - D;
	input_int[0] = 0;
	du[0] = 0;

	real u_s, d;
	u_s = (d0u0[1] - p_min) / (p_max - p_min) * 100;
	u_s = (u_s < 0) ? 0 : u_s;
	u_s = (u_s > 100) ? 100 : u_s;
	u_s = u_s - D;

	d = d0u0[0];

	x[0] = (x0 - x_minP) / (x_maxP - x_minP) * ((100 - D) - (S - J) / 2);

	int seed = std::chrono::system_clock::now().time_since_epoch().count() + 11311;
	std::default_random_engine gen_normal(seed);
	std::normal_distribution<real> randn(0.0, std_noise_controller);

	for (size_t i = 1; i < u.size(); i++) {
		if (simulation_type == ol) {
			input[i] = (u[i] - p_min) / (p_max - p_min) * 100;
			if (input[i] < 0)
				input[i] = 0;
			if (input[i] > 100)
				input[i] = 100;
		}
		else if (simulation_type == cl) {
			OP[i] = controller.pid(u[i], ((x[i - 1] * (x_max - x_min) / ((100 - D) - (S - J) / 2) + x_minP) - x_minP) / (x_maxP - x_minP) * 100 + randn(gen_normal), i);
			if (OP[i] > 100)
				OP[i] = 100;
			else if (OP[i] < 0)
				OP[i] = 0;
			input_int[i] = (Ts * OP[i] + tau_ip * input_int[i - 1]) / (Ts + tau_ip);
			input[i] = input_int[i];
		}
		else if (simulation_type == h_cl) {
			real SP = hydraulic_model(u[i], x[i - 1] * (x_maxP - x_minP) / ((100 - D) - (S - J) / 2) + x_minP, i);
			noise[i - 1] = randn(gen_normal);
			// Stem position controller
			OP[i] = controller.pid(SP, x[i - 1] + noise[i - 1], i);
			// IP model
			if (OP[i] > 100)
				OP[i] = 100;
			else if (OP[i] < 0)
				OP[i] = 0;
			input_int[i] = (Ts * OP[i] + tau_ip * input_int[i - 1]) / (Ts + tau_ip);
			input[i] = input_int[i];
			simulation_results.SP[i] = SP;
		}

		input[i] = input[i] - D;
		du[i] = input[i] - input[i - 1];

		if (du[i] * du[i - 1] <= 0.0 && stp == 0) {
			u_s = input[i - 1];
			stp = 1;
		}

		if (stp == 0) {
			x[i] = input[i] - d / 2 * (S - J);
			stp = 0;
		}
		else if ((-d * (input[i] - u_s)) > S) {
			d = -d;
			x[i] = input[i] - d / 2 * (S - J);
			stp = 0;
		}
		else if ((d * (input[i] - u_s)) > J) {
			x[i] = input[i] - d / 2 * (S - J);
			stp = 0;
		}
		else {
			x[i] = x[i - 1];
		}
	}

	for (size_t i = 0; i < u.size(); i++) {
		simulation_results.t[i] = i * Ts;


		simulation_results.x[i] = x[i] * (x_maxP - x_minP) / ((100 - D) - (S - J) / 2) + x_minP;
		if (simulation_type == h_cl)
			simulation_results.x[i] = simulation_results.x[i] + noise[i] / 100 * (x_max - x_min);

		if (simulation_type == ol)
			simulation_results.P[i] = u[i];
		else if (simulation_type == cl) {
			simulation_results.P[i] = input_int[i] * (p_max - p_min) / 100 + p_min;
			simulation_results.SP[i] = u[i];
		}
		else if (simulation_type == h_cl) {
			simulation_results.P[i] = input_int[i] * (p_max - p_min) / 100 + p_min;
		}
		simulation_results.OP[i] = OP[i];

		if (simulation_results.x[i] < x_minP) {
			simulation_results.x[i] = x_minP;
		}
		else if (simulation_results.x[i] > x_maxP) {
			simulation_results.x[i] = x_maxP;
		}
	}

	delete[] du;
	delete[] OP;
	delete[] x;
	delete[] input;
	delete[] input_int;
	delete[] noise;
}


template<typename real>
void ValveModel<real>::sim_he() {

	real* x = new real[u.size()]{ 0 };
	real* input = new real[u.size()]{ 0 };
	real* input_int = new real[u.size()]{ 0 };
	real* noise = new real[u.size()]{ 0 };
	real* OP = new real[u.size()]{ 0 };
	real x0 = pos0[0];
	real cum_u, u_r;

	allocate_sim_data(u.size());

	input[0] = (u[0] - p_min) / (p_max - p_min) * 100;
	input[0] = (input[0] < 0) ? 0 : input[0];
	input[0] = (input[0] > 100) ? 100 : input[0];
	input[0] = input[0] - D;
	input_int[0] = 0;


	x[0] = (x0 - x_minP) / (x_maxP - x_minP) * (100 - D);
	if (pos0[1] == 0) {
		cum_u = 0;
		u_r = 0;
	}
	else if (pos0[1] > 0) {
		cum_u = F_c;
		u_r = F_c;
	}
	else {
		cum_u = -F_c;
		u_r = -F_c;
	}

	int seed = std::chrono::system_clock::now().time_since_epoch().count() + 11311;
	std::default_random_engine gen_normal(seed);
	std::normal_distribution<real> randn(0.0, std_noise_controller);

	for (size_t i = 1; i < u.size(); i++) {
		if (simulation_type == ol) {
			input[i] = (u[i] - p_min) / (p_max - p_min) * 100;
			if (input[i] < 0)
				input[i] = 0;
			if (input[i] > 100)
				input[i] = 100;
		}
		else if (simulation_type == cl) {
			OP[i] = controller.pid(u[i], (x[i - 1] * (x_maxP - x_minP) / (100 - D) + x_minP) / (x_max - x_min) * 100 + randn(gen_normal), i);
			if (OP[i] > 100)
				OP[i] = 100;
			else if (OP[i] < 0)
				OP[i] = 0;
			input_int[i] = (Ts * OP[i] + tau_ip * input_int[i - 1]) / (Ts + tau_ip);
			input[i] = input_int[i];
		}
		else if (simulation_type == h_cl) {
			real SP = hydraulic_model(u[i], x[i - 1] * (x_maxP - x_minP) / (100 - D) + x_minP, i);
			noise[i - 1] = randn(gen_normal);
			// Stem position controller
			OP[i] = controller.pid(SP, x[i - 1] + noise[i - 1], i);
			// IP model
			if (OP[i] > 100)
				OP[i] = 100;
			else if (OP[i] < 0)
				OP[i] = 0;
			input_int[i] = (Ts * OP[i] + tau_ip * input_int[i - 1]) / (Ts + tau_ip);
			input[i] = input_int[i];
			simulation_results.SP[i] = SP;
		}

		input[i] = input[i] - D;

		cum_u = u_r + input[i] - input[i - 1];
		if (std::abs(cum_u) > F_s) {
			x[i] = input[i] - signal_fnc(cum_u - F_s) * F_c;
			u_r = signal_fnc(cum_u - F_s) * F_c;
		}
		else {
			x[i] = x[i - 1];
			u_r = cum_u;
		}
	}

	for (size_t i = 0; i < u.size(); i++) {
		simulation_results.t[i] = i * Ts;

		simulation_results.x[i] = x[i] * (x_maxP - x_minP) / (100 - D - F_c) + x_minP;
		if (simulation_type == h_cl)
			simulation_results.x[i] = simulation_results.x[i] + noise[i] / 100 * (x_max - x_min);

		if (simulation_type == ol)
			simulation_results.P[i] = u[i];
		else if (simulation_type == cl) {
			simulation_results.P[i] = input_int[i] * (p_max - p_min) / 100 + p_min;
			simulation_results.SP[i] = u[i];
		}
		else if (simulation_type == h_cl) {
			simulation_results.P[i] = input_int[i] * (p_max - p_min) / 100 + p_min;
		}
		simulation_results.OP[i] = OP[i];

		if (simulation_results.x[i] < x_minP) {
			simulation_results.x[i] = x_minP;
		}
		else if (simulation_results.x[i] > x_maxP) {
			simulation_results.x[i] = x_maxP;
		}
	}

	delete[] OP;
	delete[] x;
	delete[] input;
	delete[] input_int;
	delete[] noise;
}


template<typename real>
void ValveModel<real>::sim_choudhury() {

	real* v_u = new real[u.size()]{ 0 };
	real* x = new real[u.size()]{ 0 };
	real* input = new real[u.size()]{ 0 };
	real* input_int = new real[u.size()]{ 0 };
	real* noise = new real[u.size()]{ 0 };
	real* OP = new real[u.size()]{ 0 };
	real x0 = pos0[0];

	allocate_sim_data(u.size());

	int I = 0;
	input[0] = (u[0] - p_min) / (p_max - p_min) * 100;
	input[0] = (input[0] < 0) ? 0 : input[0];
	input[0] = (input[0] > 100) ? 100 : input[0];
	input[0] = input[0] - D;
	input_int[0] = 0;

	real u_s;
	u_s = (d0u0[1] - p_min) / (p_max - p_min) * 100;
	u_s = (u_s < 0) ? 0 : u_s;
	u_s = (u_s > 100) ? 100 : u_s;
	u_s = u_s - D;

	x[0] = (x0 - x_minP) / (x_maxP - x_minP) * ((100 - D) - (S - J) / 2);

	int seed = std::chrono::system_clock::now().time_since_epoch().count() + 11311;
	std::default_random_engine gen_normal(seed);
	std::normal_distribution<real> randn(0.0, std_noise_controller);

	for (size_t i = 1; i < u.size(); i++) {
		if (simulation_type == ol) {
			input[i] = (u[i] - p_min) / (p_max - p_min) * 100;
			input[i] = (input[i] < 0) ? 0 : input[i];
			input[i] = (input[i] > 100) ? 100 : input[i];
		}
		else if (simulation_type == cl) {
			OP[i] = controller.pid(u[i], ((x[i - 1] * (x_max - x_min) / ((100 - D) - (S - J) / 2) + x_minP) - x_minP) / (x_maxP - x_minP) * 100 + randn(gen_normal), i);
			if (OP[i] > 100)
				OP[i] = 100;
			else if (OP[i] < 0)
				OP[i] = 0;
			input_int[i] = (Ts * OP[i] + tau_ip * input_int[i - 1]) / (Ts + tau_ip);
			input[i] = input_int[i];
		}
		else if (simulation_type == h_cl) {
			real SP = hydraulic_model(u[i], x[i - 1] * (x_maxP - x_minP) / ((100 - D) - (S - J) / 2) + x_minP, i);
			noise[i - 1] = randn(gen_normal);
			// Stem position controller
			OP[i] = controller.pid(SP, x[i - 1] + noise[i - 1], i);
			// IP model
			if (OP[i] > 100)
				OP[i] = 100;
			else if (OP[i] < 0)
				OP[i] = 0;
			input_int[i] = (Ts * OP[i] + tau_ip * input_int[i - 1]) / (Ts + tau_ip);
			input[i] = input_int[i];
			simulation_results.SP[i] = SP;
		}

		input[i] = input[i] - D;
		v_u[i] = (input[i] - input[i - 1]) / Ts;

		if (signal_fnc(v_u[i]) == signal_fnc(v_u[i - 1])) {
			if (I == 1) {
				if (signal_fnc(input[i] - u_s) * signal_fnc(u_s - x[i - 1]) == 1.0) {
					if (std::abs(input[i] - u_s) > J) {
						I = 0;
						x[i] = input[i] - signal_fnc(v_u[i]) * (S - J) / 2;
					}
					else {
						I = 1;
						x[i] = x[i - 1];
					}
				}
				else {
					if (std::abs(input[i] - x[i - 1]) > (S + J) / 2) {
						I = 0;
						x[i] = input[i] - signal_fnc(v_u[i]) * (S - J) / 2;
					}
					else {
						I = 1;
						x[i] = x[i - 1];
					}
				}
			}
			else {
				I = 0;
				x[i] = input[i] - signal_fnc(v_u[i]) * (S - J) / 2;
			}
		}
		else {
			if (I == 0)
				u_s = input[i - 1];
			if (signal_fnc(v_u[i]) == 0.0) {
				I = 1;
				x[i] = x[i - 1];
			}
			else {
				if (std::abs(input[i] - x[i - 1]) > (S + J) / 2) {
					I = 0;
					x[i] = input[i] - signal_fnc(v_u[i]) * (S - J) / 2;
				}
				else {
					I = 1;
					x[i] = x[i - 1];
				}
			}
		}
	}

	for (size_t i = 0; i < u.size(); i++) {
		simulation_results.t[i] = i * Ts;

		simulation_results.x[i] = x[i] * (x_maxP - x_minP) / ((100 - D) - (S - J) / 2) + x_minP;
		if (simulation_type == h_cl)
			simulation_results.x[i] = simulation_results.x[i] + noise[i] / 100 * (x_max - x_min);

		if (simulation_type == ol)
			simulation_results.P[i] = u[i];
		else if (simulation_type == cl) {
			simulation_results.P[i] = input_int[i] * (p_max - p_min) / 100 + p_min;
			simulation_results.SP[i] = u[i];
		}
		else if (simulation_type == h_cl) {
			simulation_results.P[i] = input_int[i] * (p_max - p_min) / 100 + p_min;
		}
		simulation_results.OP[i] = OP[i];

		if (simulation_results.x[i] < x_minP) {
			simulation_results.x[i] = x_minP;
		}
		else if (simulation_results.x[i] > x_maxP) {
			simulation_results.x[i] = x_maxP;
		}
	}

	delete[] v_u;
	delete[] OP;
	delete[] x;
	delete[] input;
	delete[] input_int;
	delete[] noise;
}


template <typename real>
void ValveModel<real>::sim_karnopp() {

	if (Ts < dt)
		throw "Error: simulation sampling time has to be lower than data sampling time";

	size_t nSamp = size_t(Ts / dt);

	// OP initialization
	real* OP_exc = new real[nSamp + 2]{ 0.0 };
	real* P = new real[nSamp + 2]{ 0.0 };

	// Valve model with Karnopp friction model
	// valve stem position initialization
	real* x = new real[nSamp + 2]{ 0 };
	x[nSamp + 1] = pos0[0];
	x[nSamp] = pos0[0];

	// diaphragm pressure initialization
	real* P_exc = new real[nSamp + 2]{ 0 };
	real* P_us = new real[nSamp]{ 0 };

	// valve stem velocity initialization
	real* v = new real[nSamp + 2]{ 0 };

	// valve stem acceleration initialization
	real* a = new real[nSamp + 2]{ 0 };

	// friction force initialization
	real* F_at = new real[nSamp + 2]{ 0 };

	// resultant force initialization
	real* F_res = new real[nSamp + 2]{ 0 };


	allocate_sim_data(u.size());

	bool stick = true;
	real F_r{ 0 }, P_exc_old{ 0 }, P_exc_old_old{ 0 }, Fs_m_Fc{ F_s - F_c }, one_m{ 1 / m }, dt_2{ dt / 2 };

	int seed = std::chrono::system_clock::now().time_since_epoch().count() + 11311;
	std::default_random_engine gen_normal(seed);
	std::normal_distribution<real> randn(0.0, std_noise_controller);

	for (size_t i = 0; i < u.size(); i++) {

		if (i == u.size() - 1) {
			simulation_results.x[i] = x[nSamp + 1];
			if (simulation_type == ol) {
				simulation_results.P[i] = u[u.size() - 1];
			}
			else if (simulation_type == cl) {
				simulation_results.P[i] = P_exc[nSamp + 1];
				simulation_results.OP[i] = OP_exc[nSamp + 1];
				simulation_results.SP[i] = u[i];
			}
			else if (simulation_type == h_cl) {
				simulation_results.P[i] = P_exc[nSamp + 1];
				simulation_results.OP[i] = OP_exc[nSamp + 1];
				simulation_results.SP[i] = simulation_results.SP[i - 1];
			}
			simulation_results.t[i] = i * Ts + t0;
			break;
		}

		// reinitialize intermediate vectors
		if (simulation_type == ol) {
			upsample(u[i], u[i + 1], nSamp, &P_us[0]);
			P_exc_old = P_exc[nSamp + 1];
			P_exc_old_old = P_exc[nSamp];
			if (i == 0) {
				P_exc[0] = u[0];
				P_exc[1] = u[0];
			}
			else {
				P_exc[0] = P_exc_old_old;
				P_exc[1] = P_exc_old;
			}
			for (size_t ct = 0; ct < nSamp; ct++) {
				P_exc[ct + 2] = P_us[ct];
			}
		}
		else if (simulation_type == cl) {
			real OP{ 0 };
			if (i == 0) {
				OP_exc[0] = 0;
				OP_exc[1] = 0;
			}
			else {
				OP = controller.pid(u[i], (simulation_results.x[i - 1] - x_min) / (x_max - x_min) * 100 + randn(gen_normal), i);
				if (OP > 100)
					OP = 100;
				else if (OP < 0)
					OP = 0;
				OP_exc[0] = simulation_results.OP[i - 1];
				OP_exc[1] = simulation_results.OP[i - 1];
			}
			for (size_t ct = 2; ct < nSamp + 2; ++ct) {
				OP_exc[ct] = OP;
				P[ct] = (dt * (p_max - p_min) / 100 * OP_exc[ct - 1] + tau_ip * P[ct - 1]) / (dt + tau_ip);
				P_exc[ct] = P[ct] + p_min;
			}
			P[0] = P[nSamp];
			P[1] = P[nSamp + 1];
			P_exc[0] = P_exc[nSamp];
			P_exc[1] = P_exc[nSamp + 1];
			simulation_results.OP[i] = OP;
			simulation_results.SP[i] = u[i];
		}
		else if (simulation_type == h_cl) {
			real SP{ 0.0 };
			if (i > 0) {
				SP = hydraulic_model(u[i], simulation_results.x[i - 1], i);
				simulation_results.x[i - 1] += randn(gen_normal) / 100 * (x_max - x_min);
			}
			else
				SP = 0;

			// Stem position controller
			real OP{ 0 };
			if (i == 0) {
				OP_exc[0] = 0;
				OP_exc[1] = 0;
			}
			else {
				OP = controller.pid(SP, (simulation_results.x[i - 1] - x_min) / (x_max - x_min) * 100, i);
				if (OP > 100)
					OP = 100;
				else if (OP < 0)
					OP = 0;
				OP_exc[0] = simulation_results.OP[i - 1];
				OP_exc[1] = simulation_results.OP[i - 1];
			}
			for (size_t ct = 2; ct < nSamp + 2; ++ct) {
				OP_exc[ct] = OP;
				P[ct] = (dt * (p_max - p_min) / 100 * OP_exc[ct - 1] + tau_ip * P[ct - 1]) / (dt + tau_ip);
				P_exc[ct] = P[ct] + p_min;
			}
			P[0] = P[nSamp];
			P[1] = P[nSamp + 1];
			P_exc[0] = P_exc[nSamp];
			P_exc[1] = P_exc[nSamp + 1];
			simulation_results.OP[i] = OP;
			simulation_results.SP[i] = SP;
		}

		x[0] = x[nSamp];
		x[1] = x[nSamp + 1];
		v[0] = v[nSamp];
		v[1] = v[nSamp + 1];
		a[0] = a[nSamp];
		a[1] = a[nSamp + 1];
		F_at[0] = F_at[nSamp];
		F_at[1] = F_at[nSamp + 1];
		F_res[0] = F_res[nSamp];
		F_res[1] = F_res[nSamp + 1];

		simulation_results.x[i] = x[1];
		simulation_results.P[i] = P_exc[1];
		simulation_results.t[i] = real(i * Ts) + t0;

		if (isnan(simulation_results.x[i])) {
			delete[] x;
			delete[] P_exc;
			delete[] P_us;
			delete[] v;
			delete[] a;
			delete[] F_at;
			delete[] F_res;
			delete[] OP_exc;
			delete[] P;
			clear_sim_data();
			return;
		}

		for (size_t j = 2; j < nSamp + 2; j++) {
			if (v[j - 2] == 0 && !stick)
				stick = false;
			else if (v[j - 2] * v[j - 1] <= 0)
				stick = true;
			else
				stick = false;
			if (stick) {
				F_r = S_a * P_exc[j - 1] - k * x[j - 1] - F_init;
				real absFr = (F_r < 0 ? -F_r : F_r);
				if (F_r > 0)
					F_at[j] = (absFr < F_s ? absFr : F_s);
				else if (F_r < 0)
					F_at[j] = (absFr < F_s ? -absFr : -F_s);
				else
					F_at[j] = 0;

				F_res[j] = F_r - F_at[j];
				v[j - 1] = 0;

				real absFat = (F_at[j] < 0 ? -F_at[j] : F_at[j]);
				stick = (absFat >= F_s ? false : true);
			}
			else {
				if (v[j - 1] > 0) {
					real v_vs = v[j - 1] / v_s;
					real exp_v2 = std::exp(-v_vs * v_vs);
					F_at[j] = (F_c + Fs_m_Fc * exp_v2) + F_v * v[j - 1];
				}
				else if (v[j - 1] < 0) {
					real v_vs = v[j - 1] / v_s;
					real exp_v2 = std::exp(-v_vs * v_vs);
					F_at[j] = -(F_c + Fs_m_Fc * exp_v2) + F_v * v[j - 1];
				}
				else {
					F_at[j] = 0;
				}
				F_res[j] = S_a * P_exc[j - 1] - k * x[j - 1] - F_at[j] - F_init;
			}

			a[j] = one_m * F_res[j];
			v[j] = dt_2 * (a[j] + a[j - 1]) + v[j - 1];
			x[j] = dt_2 * (v[j] + v[j - 1]) + x[j - 1];

			if (x[j] < x_min) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_min;
			}
			else if (x[j] > x_max) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_max;
			}
		}

	}
	delete[] x;
	delete[] P_exc;
	delete[] P_us;
	delete[] v;
	delete[] a;
	delete[] F_at;
	delete[] F_res;
	delete[] OP_exc;
	delete[] P;
}


template<typename real>
void ValveModel<real>::sim_lugre() {
	if (Ts < dt)
		throw "Error: simulation sampling time has to be lower than data sampling time";

	size_t nSamp = size_t(Ts / dt);

	// OP initialization
	real* OP_exc = new real[nSamp + 2]{ 0.0 };
	real* P = new real[nSamp + 2]{ 0.0 };

	// Valve model with Karnopp friction model
	// valve stem position initialization
	real* x = new real[nSamp + 2]{ 0 };
	x[nSamp + 1] = pos0[0];
	x[nSamp] = pos0[0];

	// diaphragm pressure initialization
	real* P_exc = new real[nSamp + 2]{ 0 };
	real* P_us = new real[nSamp]{ 0 };

	// valve stem velocity initialization
	real* v = new real[nSamp + 2]{ 0 };

	// valve stem acceleration initialization
	real* a = new real[nSamp + 2]{ 0 };

	// friction force initialization
	real* F_at = new real[nSamp + 2]{ 0 };

	// resultant force initialization
	real* F_res = new real[nSamp + 2]{ 0 };

	// internal states vector
	real* z = new real[nSamp + 2]{ 0 };
	real z0;

	if (pos0[1] == 1)
		z0 = std::abs(S_a * u[0] - k * pos0[0] - F_init) / sigma_0;
	else if (pos0[1] == -1)
		z0 = -std::abs(S_a * u[0] - k * pos0[0] - F_init) / sigma_0;
	else
		z0 = 0.0;
	z[nSamp + 1] = z0;
	z[nSamp] = z0;

	allocate_sim_data(u.size());

	real P_exc_old{ 0 }, P_exc_old_old{ 0 };

	real g_v{ 0 }, dot_z{ 0 }, dot_z_ant{ 0 }, Fs_m_Fc{ F_s - F_c }, dt_u_2{ dt / 2 }, inv_m{ 1 / m }, inv_sigma0{ 1 / sigma_0 };

	int seed = std::chrono::system_clock::now().time_since_epoch().count() + 11311;
	std::default_random_engine gen_normal(seed);
	std::normal_distribution<real> randn(0.0, std_noise_controller);

	for (size_t i = 0; i < u.size(); i++) {

		if (i == u.size() - 1) {
			simulation_results.x[i] = x[nSamp + 1];
			if (simulation_type == ol) {
				simulation_results.P[i] = u[u.size() - 1];
			}
			else if (simulation_type == cl) {
				simulation_results.P[i] = P_exc[nSamp + 1];
				simulation_results.OP[i] = OP_exc[nSamp + 1];
				simulation_results.SP[i] = u[i];
			}
			else if (simulation_type == h_cl) {
				simulation_results.P[i] = P_exc[nSamp + 1];
				simulation_results.OP[i] = OP_exc[nSamp + 1];
				simulation_results.SP[i] = simulation_results.SP[i - 1];
			}
			simulation_results.t[i] = i * Ts + t0;
			break;
		}

		// reinitialize intermediate vectors
		if (simulation_type == ol) {
			upsample(u[i], u[i + 1], nSamp, &P_us[0]);
			P_exc_old = P_exc[nSamp + 1];
			P_exc_old_old = P_exc[nSamp];
			if (i == 0) {
				P_exc[0] = u[0];
				P_exc[1] = u[0];
			}
			else {
				P_exc[0] = P_exc_old_old;
				P_exc[1] = P_exc_old;
			}
			for (size_t ct = 0; ct < nSamp; ct++) {
				P_exc[ct + 2] = P_us[ct];
			}
		}
		else if (simulation_type == cl) {
			real OP{ 0 };
			if (i == 0) {
				OP_exc[0] = 0;
				OP_exc[1] = 0;
			}
			else {
				OP = controller.pid(u[i], (simulation_results.x[i - 1] - x_min) / (x_max - x_min) * 100 + randn(gen_normal), i);
				if (OP > 100)
					OP = 100;
				else if (OP < 0)
					OP = 0;
				OP_exc[0] = simulation_results.OP[i - 1];
				OP_exc[1] = simulation_results.OP[i - 1];
			}
			for (size_t ct = 2; ct < nSamp + 2; ++ct) {
				OP_exc[ct] = OP;
				P[ct] = (dt * (p_max - p_min) / 100 * OP_exc[ct - 1] + tau_ip * P[ct - 1]) / (dt + tau_ip);
				P_exc[ct] = P[ct] + p_min;
			}
			P[0] = P[nSamp];
			P[1] = P[nSamp + 1];
			P_exc[0] = P_exc[nSamp];
			P_exc[1] = P_exc[nSamp + 1];
			simulation_results.OP[i] = OP;
			simulation_results.SP[i] = u[i];
		}
		else if (simulation_type == h_cl) {
			real SP{ 0.0 };
			if (i > 0) {
				SP = hydraulic_model(u[i], simulation_results.x[i - 1], i);
				simulation_results.x[i - 1] += randn(gen_normal) / 100 * (x_max - x_min);
			}
			else
				SP = 0;

			// Stem position controller
			real OP{ 0 };
			if (i == 0) {
				OP_exc[0] = 0;
				OP_exc[1] = 0;
			}
			else {
				OP = controller.pid(SP, (simulation_results.x[i - 1] - x_min) / (x_max - x_min) * 100, i);
				if (OP > 100)
					OP = 100;
				else if (OP < 0)
					OP = 0;
				OP_exc[0] = simulation_results.OP[i - 1];
				OP_exc[1] = simulation_results.OP[i - 1];
			}
			for (size_t ct = 2; ct < nSamp + 2; ++ct) {
				OP_exc[ct] = OP;
				P[ct] = (dt * (p_max - p_min) / 100 * OP_exc[ct - 1] + tau_ip * P[ct - 1]) / (dt + tau_ip);
				P_exc[ct] = P[ct] + p_min;
			}
			P[0] = P[nSamp];
			P[1] = P[nSamp + 1];
			P_exc[0] = P_exc[nSamp];
			P_exc[1] = P_exc[nSamp + 1];
			simulation_results.OP[i] = OP;
			simulation_results.SP[i] = SP;
		}

		x[0] = x[nSamp];
		x[1] = x[nSamp + 1];
		v[0] = v[nSamp];
		v[1] = v[nSamp + 1];
		a[0] = a[nSamp];
		a[1] = a[nSamp + 1];
		F_at[0] = F_at[nSamp];
		F_at[1] = F_at[nSamp + 1];
		F_res[0] = F_res[nSamp];
		F_res[1] = F_res[nSamp + 1];
		z[0] = z[nSamp];
		z[1] = z[nSamp + 1];

		simulation_results.x[i] = x[1];
		simulation_results.P[i] = P_exc[1];
		simulation_results.t[i] = i * Ts + t0;

		if (isnan(simulation_results.x[i])) {
			delete[] x;
			delete[] P_exc;
			delete[] P_us;
			delete[] v;
			delete[] a;
			delete[] F_at;
			delete[] F_res;
			delete[] z;
			delete[] OP_exc;
			delete[] P;
			clear_sim_data();
			return;
		}

		for (size_t j = 2; j < nSamp + 2; j++) {
			real v_vs = v[j - 1] / v_s;
			real exp_vs = std::exp(-v_vs * v_vs);
			g_v = inv_sigma0 * (F_c + Fs_m_Fc * exp_vs);

			real abs_v = (v[j - 1] < 0 ? -v[j - 1] : v[j - 1]);
			dot_z = v[j - 1] - abs_v / g_v * z[j - 1];

			z[j] = dt_u_2 * (dot_z + dot_z_ant) + z[j - 1];

			F_at[j] = sigma_0 * z[j] + sigma_1 * dot_z + F_v * v[j - 1];

			F_res[j] = S_a * P_exc[j - 1] - k * x[j - 1] - F_init - F_at[j];

			a[j] = inv_m * F_res[j];
			v[j] = dt_u_2 * (a[j] + a[j - 1]) + v[j - 1];
			x[j] = dt_u_2 * (v[j] + v[j - 1]) + x[j - 1];

			if (x[j] < x_min) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_min;
			}
			else if (x[j] > x_max) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_max;
			}

			dot_z_ant = dot_z;
		}

	}
	delete[] x;
	delete[] P_exc;
	delete[] P_us;
	delete[] v;
	delete[] a;
	delete[] F_at;
	delete[] F_res;
	delete[] z;
	delete[] OP_exc;
	delete[] P;
}


template<typename real>
void ValveModel<real>::sim_gms() {
	if (Ts < dt)
		throw "Error: simulation sampling time has to be lower than data sampling time";

	size_t nSamp = size_t(Ts / dt);

	// OP initialization
	real* OP_exc = new real[nSamp + 2]{ 0.0 };
	real* P = new real[nSamp + 2]{ 0.0 };

	// Valve model with Karnopp friction model
	// valve stem position initialization
	real* x = new real[nSamp + 2]{ 0 };
	x[nSamp + 1] = pos0[0];
	x[nSamp] = pos0[0];

	// diaphragm pressure initialization
	real* P_exc = new real[nSamp + 2]{ 0 };
	real* P_us = new real[nSamp]{ 0 };

	// valve stem velocity initialization
	real* v = new real[nSamp + 2]{ 0 };

	// valve stem acceleration initialization
	real* a = new real[nSamp + 2]{ 0 };

	// friction force initialization
	real* F_at = new real[nSamp + 2]{ 0 };

	// resultant force initialization
	real* F_res = new real[nSamp + 2]{ 0 };

	// internal states vector
	real* z_1 = new real[nSamp + 2]{ 0 };
	real* z_2 = new real[nSamp + 2]{ 0 };
	real* z_3 = new real[nSamp + 2]{ 0 };
	real z_1_0{ 0 }, z_2_0{ 0 }, z_3_0{ 0 };


	if (pos0[1] == 1) {
		z_1_0 = std::abs(S_a * u[0] - k * pos0[0] - F_init) * alpha_1 / kappa_1;
		z_2_0 = std::abs(S_a * u[0] - k * pos0[0] - F_init) * alpha_2 / kappa_2;
		z_3_0 = std::abs(S_a * u[0] - k * pos0[0] - F_init) * alpha_3 / kappa_3;
	}
	else if (pos0[1] == -1) {
		z_1_0 = -std::abs(S_a * u[0] - k * pos0[0] - F_init) * alpha_1 / kappa_1;
		z_2_0 = -std::abs(S_a * u[0] - k * pos0[0] - F_init) * alpha_2 / kappa_2;
		z_3_0 = -std::abs(S_a * u[0] - k * pos0[0] - F_init) * alpha_3 / kappa_3;
	}
	else {
		z_1_0 = 0;
		z_2_0 = 0;
		z_3_0 = 0;
	}
	z_1[nSamp + 1] = z_1_0;
	z_1[nSamp] = z_1_0;
	z_2[nSamp + 1] = z_2_0;
	z_2[nSamp] = z_2_0;
	z_3[nSamp + 1] = z_3_0;
	z_3[nSamp] = z_3_0;

	allocate_sim_data(u.size());

	real P_exc_old{ 0 }, P_exc_old_old{ 0 };

	real F_at_1{ 0 }, F_at_1_ant{ z_1_0 * kappa_1 * alpha_1 }, F_at_2{ 0 }, F_at_2_ant{ z_2_0 * kappa_2 * alpha_2 }, F_at_3{ 0 }, F_at_3_ant{ z_3_0 * kappa_3 * alpha_3 };
	real dot_z_1{ 0 }, dot_z_1_ant{ 0 }, dot_z_2{ 0 }, dot_z_2_ant{ 0 }, dot_z_3{ 0 }, dot_z_3_ant{ 0 };
	real sinal{ 0 }, s_v{ 0 }, inv_m{ 1 / m }, dt_u_2{ dt / 2 }, Fs_m_Fc{ F_s - F_c };
	real C_alpha1_u_kappa1{ C * alpha_1 / kappa_1 }, C_alpha2_u_kappa2{ C * alpha_2 / kappa_2 }, C_alpha3_u_kappa3{ C * alpha_3 / kappa_3 };
	bool stick_1{ true }, stick_2{ true }, stick_3{ true };
	//__m256d v_kappa = _mm256_set_pd(kappa_1, kappa_2, kappa_3, 0.0f);
	//__m256d v_nu = _mm256_set_pd(nu_1, nu_2, nu_3, 0.0f);


	int seed = std::chrono::system_clock::now().time_since_epoch().count() + 11311;
	std::default_random_engine gen_normal(seed);
	std::normal_distribution<real> randn(0.0, std_noise_controller);

	for (size_t i = 0; i < u.size(); i++) {

		if (i == u.size() - 1) {
			simulation_results.x[i] = x[nSamp + 1];
			if (simulation_type == ol) {
				simulation_results.P[i] = u[u.size() - 1];
			}
			else if (simulation_type == cl) {
				simulation_results.P[i] = P_exc[nSamp + 1];
				simulation_results.OP[i] = OP_exc[nSamp + 1];
				simulation_results.SP[i] = u[i];
			}
			else if (simulation_type == h_cl) {
				simulation_results.P[i] = P_exc[nSamp + 1];
				simulation_results.OP[i] = OP_exc[nSamp + 1];
				simulation_results.SP[i] = simulation_results.SP[i - 1];
			}
			simulation_results.t[i] = i * Ts + t0;
			break;
		}


		// reinitialize intermediate vectors
		if (simulation_type == ol) {
			upsample(u[i], u[i + 1], nSamp, &P_us[0]);
			P_exc_old = P_exc[nSamp + 1];
			P_exc_old_old = P_exc[nSamp];
			if (i == 0) {
				P_exc[0] = u[0];
				P_exc[1] = u[0];
			}
			else {
				P_exc[0] = P_exc_old_old;
				P_exc[1] = P_exc_old;
			}
			for (size_t ct = 0; ct < nSamp; ct++) {
				P_exc[ct + 2] = P_us[ct];
			}
		}
		else if (simulation_type == cl) {
			real OP{ 0 };
			if (i == 0) {
				OP_exc[0] = 0;
				OP_exc[1] = 0;
			}
			else {
				OP = controller.pid(u[i], (simulation_results.x[i - 1] - x_min) / (x_max - x_min) * 100 + randn(gen_normal), i);
				if (OP > 100)
					OP = 100;
				else if (OP < 0)
					OP = 0;
				OP_exc[0] = simulation_results.OP[i - 1];
				OP_exc[1] = simulation_results.OP[i - 1];
			}
			for (size_t ct = 2; ct < nSamp + 2; ++ct) {
				OP_exc[ct] = OP;
				P[ct] = (dt * (p_max - p_min) / 100 * OP_exc[ct - 1] + tau_ip * P[ct - 1]) / (dt + tau_ip);
				P_exc[ct] = P[ct] + p_min;
			}
			P[0] = P[nSamp];
			P[1] = P[nSamp + 1];
			P_exc[0] = P_exc[nSamp];
			P_exc[1] = P_exc[nSamp + 1];
			simulation_results.OP[i] = OP;
			simulation_results.SP[i] = u[i];
		}
		else if (simulation_type == h_cl) {
			real SP{ 0.0 };
			if (i > 0) {
				SP = hydraulic_model(u[i], simulation_results.x[i - 1], i);
				simulation_results.x[i - 1] += randn(gen_normal) / 100 * (x_max - x_min);
			}
			else
				SP = 0;

			// Stem position controller
			real OP{ 0 };
			if (i == 0) {
				OP_exc[0] = 0;
				OP_exc[1] = 0;
			}
			else {
				OP = controller.pid(SP, (simulation_results.x[i - 1] - x_min) / (x_max - x_min) * 100, i);
				if (OP > 100)
					OP = 100;
				else if (OP < 0)
					OP = 0;
				OP_exc[0] = simulation_results.OP[i - 1];
				OP_exc[1] = simulation_results.OP[i - 1];
			}
			for (size_t ct = 2; ct < nSamp + 2; ++ct) {
				OP_exc[ct] = OP;
				P[ct] = (dt * (p_max - p_min) / 100 * OP_exc[ct - 1] + tau_ip * P[ct - 1]) / (dt + tau_ip);
				P_exc[ct] = P[ct] + p_min;
			}
			P[0] = P[nSamp];
			P[1] = P[nSamp + 1];
			P_exc[0] = P_exc[nSamp];
			P_exc[1] = P_exc[nSamp + 1];
			simulation_results.OP[i] = OP;
			simulation_results.SP[i] = SP;
		}


		x[0] = x[nSamp];
		x[1] = x[nSamp + 1];
		v[0] = v[nSamp];
		v[1] = v[nSamp + 1];
		a[0] = a[nSamp];
		a[1] = a[nSamp + 1];
		F_at[0] = F_at[nSamp];
		F_at[1] = F_at[nSamp + 1];
		F_res[0] = F_res[nSamp];
		F_res[1] = F_res[nSamp + 1];
		z_1[0] = z_1[nSamp];
		z_1[1] = z_1[nSamp + 1];
		z_2[0] = z_2[nSamp];
		z_2[1] = z_2[nSamp + 1];
		z_3[0] = z_3[nSamp];
		z_3[1] = z_3[nSamp + 1];

		simulation_results.x[i] = x[1];
		simulation_results.P[i] = P_exc[1];
		simulation_results.t[i] = i * Ts + t0;

		if (isnan(simulation_results.x[i])) {
			delete[] x;
			delete[] P_exc;
			delete[] P_us;
			delete[] v;
			delete[] a;
			delete[] F_at;
			delete[] F_res;
			delete[] z_1;
			delete[] z_2;
			delete[] z_3;
			delete[] OP_exc;
			delete[] P;
			clear_sim_data();
			return;
		}


		for (size_t j = 2; j < nSamp + 2; j++) {
			if (v[j - 1] > 0)
				sinal = 1;
			else if (v[j - 1] < 0)
				sinal = -1;
			else
				sinal = 0;

			if (v[j - 2] == 0.0 && stick_1 == 0)
				stick_1 = false;
			else if (v[j - 2] * v[j - 1] <= 0)
				stick_1 = true;
			if (v[j - 2] == 0.0 && stick_2 == 0)
				stick_2 = false;
			else if (v[j - 2] * v[j - 1] <= 0)
				stick_2 = true;
			if (v[j - 2] == 0.0 && stick_3 == 0)
				stick_3 = false;
			else if (v[j - 2] * v[j - 1] <= 0)
				stick_3 = true;

			real v_vs = v[j - 1] / v_s;
			real exp_vs = std::exp(-v_vs * v_vs);
			s_v = (F_c + Fs_m_Fc * exp_vs);

			// First element
			if (stick_1) {
				dot_z_1 = v[j - 1];
				z_1[j] = dt_u_2 * (v[j - 1] + v[j - 2]) + z_1[j - 1];
			}
			else {
				dot_z_1 = C_alpha1_u_kappa1 * (sinal - F_at_1_ant / alpha_1 / s_v);
				z_1[j] = dt_u_2 * (dot_z_1 + dot_z_1_ant) + z_1[j - 1];
			}

			// Second element
			if (stick_2) {
				dot_z_2 = v[j - 1];
				z_2[j] = dt_u_2 * (v[j - 1] + v[j - 2]) + z_2[j - 1];
			}
			else {
				dot_z_2 = C_alpha2_u_kappa2 * (sinal - F_at_2_ant / alpha_2 / s_v);
				z_2[j] = dt_u_2 * (dot_z_2 + dot_z_2_ant) + z_2[j - 1];
			}

			// Third element
			if (stick_3) {
				dot_z_3 = v[j - 1];
				z_3[j] = dt_u_2 * (v[j - 1] + v[j - 2]) + z_3[j - 1];
			}
			else {
				dot_z_3 = C_alpha3_u_kappa3 * (sinal - F_at_3_ant / alpha_3 / s_v);
				z_3[j] = dt_u_2 * (dot_z_3 + dot_z_3_ant) + z_3[j - 1];
			}

			F_at_1 = kappa_1 * z_1[j] + nu_1 * dot_z_1;
			F_at_2 = kappa_2 * z_2[j] + nu_2 * dot_z_2;
			F_at_3 = kappa_3 * z_3[j] + nu_3 * dot_z_3;

			real absFat1 = (F_at_1 < 0 ? -F_at_1 : F_at_1);
			real absFat2 = (F_at_2 < 0 ? -F_at_2 : F_at_2);
			real absFat3 = (F_at_3 < 0 ? -F_at_3 : F_at_3);
			if (absFat1 > alpha_1 * s_v)
				stick_1 = false;
			if (absFat2 > alpha_2 * s_v)
				stick_2 = false;
			if (absFat3 > alpha_3 * s_v)
				stick_3 = false;

			F_at[j] = F_at_1 + F_at_2 + F_at_3 + F_v * v[j - 1];

			F_res[j] = S_a * P_exc[j - 1] - k * x[j - 1] - F_init - F_at[j];

			a[j] = inv_m * F_res[j];
			v[j] = dt_u_2 * (a[j] + a[j - 1]) + v[j - 1];
			x[j] = dt_u_2 * (v[j] + v[j - 1]) + x[j - 1];

			if (x[j] < x_min) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_min;
			}
			else if (x[j] > x_max) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_max;
			}

			F_at_1_ant = F_at_1;
			F_at_2_ant = F_at_2;
			F_at_3_ant = F_at_3;
			dot_z_1_ant = dot_z_1;
			dot_z_2_ant = dot_z_2;
			dot_z_3_ant = dot_z_3;
		}

	}
	delete[] x;
	delete[] P_exc;
	delete[] P_us;
	delete[] v;
	delete[] a;
	delete[] F_at;
	delete[] F_res;
	delete[] z_1;
	delete[] z_2;
	delete[] z_3;
	delete[] OP_exc;
	delete[] P;
}


template<typename real>
void ValveModel<real>::sim_sgms() {
	if (Ts < dt)
		throw "Error: simulation sampling time has to be lower than data sampling time";

	size_t nSamp = size_t(Ts / dt);

	// OP initialization
	real* OP_exc = new real[nSamp + 2]{ 0.0 };
	real* P = new real[nSamp + 2]{ 0.0 };

	// Valve model with Karnopp friction model
	// valve stem position initialization
	real* x = new real[nSamp + 2]{ 0 };
	x[nSamp + 1] = pos0[0];
	x[nSamp] = pos0[0];

	// diaphragm pressure initialization
	real* P_exc = new real[nSamp + 2]{ 0 };
	real* P_us = new real[nSamp]{ 0 };

	// valve stem velocity initialization
	real* v = new real[nSamp + 2]{ 0 };

	// valve stem acceleration initialization
	real* a = new real[nSamp + 2]{ 0 };

	// friction force initialization
	real* F_at = new real[nSamp + 2]{ 0 };

	// resultant force initialization
	real* F_res = new real[nSamp + 2]{ 0 };

	// internal states vector
	real* F_at_1 = new real[nSamp + 2]{ 0 };
	real* F_at_2 = new real[nSamp + 2]{ 0 };

	real F_at_1_0, F_at_2_0;
	if (pos0[1] == 1) {
		F_at_1_0 = std::abs(S_a * u[0] - k * pos0[0] - F_init) * alpha_1;
		F_at_2_0 = std::abs(S_a * u[0] - k * pos0[0] - F_init) * alpha_2;
	}
	else if (pos0[1] == -1) {
		F_at_1_0 = -std::abs(S_a * u[0] - k * pos0[0] - F_init) * alpha_1;
		F_at_2_0 = -std::abs(S_a * u[0] - k * pos0[0] - F_init) * alpha_2;
	}
	else {
		F_at_1_0 = 0;
		F_at_2_0 = 0;
	}
	F_at_1[nSamp + 1] = F_at_1_0;
	F_at_1[nSamp] = F_at_1_0;
	F_at_2[nSamp + 1] = F_at_2_0;
	F_at_2[nSamp] = F_at_2_0;

	allocate_sim_data(u.size());

	real P_exc_old{ 0 }, P_exc_old_old{ 0 };

	real dot_F_at_1{ 0 }, dot_F_at_1_ant{ 0 }, dot_F_at_2{ 0 }, dot_F_at_2_ant{ 0 };
	real sinal{ 0 }, s_v{ 0 };
	bool stick_1{ true }, stick_2{ true };
	real inv_m{ 1 / m }, dt_u_2{ dt / 2 }, Fs_m_Fc{ F_s - F_c };


	int seed = std::chrono::system_clock::now().time_since_epoch().count() + 11311;
	std::default_random_engine gen_normal(seed);
	std::normal_distribution<real> randn(0.0, std_noise_controller);

	for (size_t i = 0; i < u.size(); i++) {

		if (i == u.size() - 1) {
			simulation_results.x[i] = x[nSamp + 1];
			if (simulation_type == ol) {
				simulation_results.P[i] = u[u.size() - 1];
			}
			else if (simulation_type == cl) {
				simulation_results.P[i] = P_exc[nSamp + 1];
				simulation_results.OP[i] = OP_exc[nSamp + 1];
				simulation_results.SP[i] = u[i];
			}
			else if (simulation_type == h_cl) {
				simulation_results.P[i] = P_exc[nSamp + 1];
				simulation_results.OP[i] = OP_exc[nSamp + 1];
				simulation_results.SP[i] = simulation_results.SP[i - 1];
			}
			simulation_results.t[i] = i * Ts + t0;
			break;
		}


		// reinitialize intermediate vectors
		if (simulation_type == ol) {
			upsample(u[i], u[i + 1], nSamp, &P_us[0]);
			P_exc_old = P_exc[nSamp + 1];
			P_exc_old_old = P_exc[nSamp];
			if (i == 0) {
				P_exc[0] = u[0];
				P_exc[1] = u[0];
			}
			else {
				P_exc[0] = P_exc_old_old;
				P_exc[1] = P_exc_old;
			}
			for (size_t ct = 0; ct < nSamp; ct++) {
				P_exc[ct + 2] = P_us[ct];
			}
		}
		else if (simulation_type == cl) {
			real OP{ 0 };
			if (i == 0) {
				OP_exc[0] = 0;
				OP_exc[1] = 0;
			}
			else {
				OP = controller.pid(u[i], (simulation_results.x[i - 1] - x_min) / (x_max - x_min) * 100 + randn(gen_normal), i);
				if (OP > 100)
					OP = 100;
				else if (OP < 0)
					OP = 0;
				OP_exc[0] = simulation_results.OP[i - 1];
				OP_exc[1] = simulation_results.OP[i - 1];
			}
			for (size_t ct = 2; ct < nSamp + 2; ++ct) {
				OP_exc[ct] = OP;
				P[ct] = (dt * (p_max - p_min) / 100 * OP_exc[ct - 1] + tau_ip * P[ct - 1]) / (dt + tau_ip);
				P_exc[ct] = P[ct] + p_min;
			}
			P[0] = P[nSamp];
			P[1] = P[nSamp + 1];
			P_exc[0] = P_exc[nSamp];
			P_exc[1] = P_exc[nSamp + 1];
			simulation_results.OP[i] = OP;
			simulation_results.SP[i] = u[i];
		}
		else if (simulation_type == h_cl) {
			real SP{ 0.0 };
			if (i > 0) {
				SP = hydraulic_model(u[i], simulation_results.x[i - 1], i);
				simulation_results.x[i - 1] += randn(gen_normal) / 100 * (x_max - x_min);
			}
			else
				SP = 0;

			// Stem position controller
			real OP{ 0 };
			if (i == 0) {
				OP_exc[0] = 0;
				OP_exc[1] = 0;
			}
			else {
				OP = controller.pid(SP, (simulation_results.x[i - 1] - x_min) / (x_max - x_min) * 100, i);
				if (OP > 100)
					OP = 100;
				else if (OP < 0)
					OP = 0;
				OP_exc[0] = simulation_results.OP[i - 1];
				OP_exc[1] = simulation_results.OP[i - 1];
			}
			for (size_t ct = 2; ct < nSamp + 2; ++ct) {
				OP_exc[ct] = OP;
				P[ct] = (dt * (p_max - p_min) / 100 * OP_exc[ct - 1] + tau_ip * P[ct - 1]) / (dt + tau_ip);
				P_exc[ct] = P[ct] + p_min;
			}
			P[0] = P[nSamp];
			P[1] = P[nSamp + 1];
			P_exc[0] = P_exc[nSamp];
			P_exc[1] = P_exc[nSamp + 1];
			simulation_results.OP[i] = OP;
			simulation_results.SP[i] = SP;
		}


		x[0] = x[nSamp];
		x[1] = x[nSamp + 1];
		v[0] = v[nSamp];
		v[1] = v[nSamp + 1];
		a[0] = a[nSamp];
		a[1] = a[nSamp + 1];
		F_at[0] = F_at[nSamp];
		F_at[1] = F_at[nSamp + 1];
		F_res[0] = F_res[nSamp];
		F_res[1] = F_res[nSamp + 1];
		F_at_1[0] = F_at_1[nSamp];
		F_at_1[1] = F_at_1[nSamp + 1];
		F_at_2[0] = F_at_2[nSamp];
		F_at_2[1] = F_at_2[nSamp + 1];

		simulation_results.x[i] = x[1];
		simulation_results.P[i] = P_exc[1];
		simulation_results.t[i] = i * Ts + t0;

		if (isnan(simulation_results.x[i])) {
			delete[] x;
			delete[] P_exc;
			delete[] P_us;
			delete[] v;
			delete[] a;
			delete[] F_at;
			delete[] F_res;
			delete[] F_at_1;
			delete[] F_at_2;
			delete[] OP_exc;
			delete[] P;
			clear_sim_data();
			return;
		}

		for (size_t j = 2; j < nSamp + 2; j++) {
			if (v[j - 1] > 0)
				sinal = 1;
			else if (v[j - 1] < 0)
				sinal = -1;
			else
				sinal = 0;

			if (v[j - 2] == 0.0 && stick_1 == 0)
				stick_1 = false;
			else if (v[j - 2] * v[j - 1] <= 0)
				stick_1 = true;

			if (v[j - 2] == 0.0 && stick_2 == 0)
				stick_2 = false;
			else if (v[j - 2] * v[j - 1] <= 0)
				stick_2 = true;

			real v_vs = v[j - 1] / v_s;
			real exp_vs = std::exp(-v_vs * v_vs);
			s_v = (F_c + Fs_m_Fc * exp_vs);

			// First element
			if (stick_1)
				dot_F_at_1 = kappa_1 * v[j - 1];
			else
				dot_F_at_1 = C * (sinal - F_at_1[j - 1] / (alpha_1 * s_v));

			// Second element
			if (stick_2)
				dot_F_at_2 = kappa_2 * v[j - 1];
			else
				dot_F_at_2 = C * (sinal - F_at_2[j - 1] / (alpha_2 * s_v));


			F_at_1[j] = dt_u_2 * (dot_F_at_1_ant + dot_F_at_1) + F_at_1[j - 1];
			F_at_2[j] = dt_u_2 * (dot_F_at_2_ant + dot_F_at_2) + F_at_2[j - 1];

			real absFat1 = (F_at_1[j] < 0 ? -F_at_1[j] : F_at_1[j]);
			real absFat2 = (F_at_2[j] < 0 ? -F_at_2[j] : F_at_2[j]);
			if (absFat1 > alpha_1 * s_v)
				stick_1 = false;
			if (absFat2 > alpha_2 * s_v)
				stick_2 = false;

			F_at[j] = F_at_1[j] + F_at_2[j] + F_v * v[j - 1];

			F_res[j] = S_a * P_exc[j - 1] - k * x[j - 1] - F_init - F_at[j];

			a[j] = inv_m * F_res[j];
			v[j] = dt_u_2 * (a[j] + a[j - 1]) + v[j - 1];
			x[j] = dt_u_2 * (v[j] + v[j - 1]) + x[j - 1];

			if (x[j] < x_min) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_min;
			}
			else if (x[j] > x_max) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_max;
			}
			dot_F_at_1_ant = dot_F_at_1;
			dot_F_at_2_ant = dot_F_at_2;
		}

	}
	delete[] x;
	delete[] P_exc;
	delete[] P_us;
	delete[] v;
	delete[] a;
	delete[] F_at;
	delete[] F_res;
	delete[] F_at_1;
	delete[] F_at_2;
	delete[] OP_exc;
	delete[] P;
}


template<typename real>
void ValveModel<real>::sim_gms1() {
	if (Ts < dt)
		throw "Error: simulation sampling time has to be lower than data sampling time";

	size_t nSamp = size_t(Ts / dt);

	// OP initialization
	real* OP_exc = new real[nSamp + 2]{ 0.0 };
	real* P = new real[nSamp + 2]{ 0.0 };

	// Valve model with Karnopp friction model
	// valve stem position initialization
	real* x = new real[nSamp + 2]{ 0 };
	x[nSamp + 1] = pos0[0];
	x[nSamp] = pos0[0];

	// diaphragm pressure initialization
	real* P_exc = new real[nSamp + 2]{ 0 };
	real* P_us = new real[nSamp]{ 0 };

	// valve stem velocity initialization
	real* v = new real[nSamp + 2]{ 0 };

	// valve stem acceleration initialization
	real* a = new real[nSamp + 2]{ 0 };

	// friction force initialization
	real* F_at = new real[nSamp + 2]{ 0 };

	// resultant force initialization
	real* F_res = new real[nSamp + 2]{ 0 };

	// internal states vector
	real* z_1 = new real[nSamp + 2]{ 0 };
	real z_1_0{ 0 };


	if (pos0[1] == 1) {
		z_1_0 = std::abs(S_a * u[0] - k * pos0[0] - F_init) / kappa_1;
	}
	else if (pos0[1] == -1) {
		z_1_0 = -std::abs(S_a * u[0] - k * pos0[0] - F_init) / kappa_1;
	}
	else {
		z_1_0 = 0;
	}
	z_1[nSamp + 1] = z_1_0;
	z_1[nSamp] = z_1_0;

	allocate_sim_data(u.size());

	real P_exc_old{ 0 }, P_exc_old_old{ 0 };

	real F_at_1{ 0 }, F_at_1_ant{ z_1_0 * kappa_1 };
	real dot_z_1{ 0 }, dot_z_1_ant{ 0 };
	real sinal{ 0 }, s_v{ 0 };
	real inv_m{ 1 / m }, dt_u_2{ dt / 2 }, Fs_m_Fc{ F_s - F_c }, C_u_kappa1{ C / kappa_1 };
	bool stick_1{ true };


	int seed = std::chrono::system_clock::now().time_since_epoch().count() + 11311;
	std::default_random_engine gen_normal(seed);
	std::normal_distribution<real> randn(0.0, std_noise_controller);

	for (size_t i = 0; i < u.size(); i++) {

		if (i == u.size() - 1) {
			simulation_results.x[i] = x[nSamp + 1];
			if (simulation_type == ol) {
				simulation_results.P[i] = u[u.size() - 1];
			}
			else if (simulation_type == cl) {
				simulation_results.P[i] = P_exc[nSamp + 1];
				simulation_results.OP[i] = OP_exc[nSamp + 1];
				simulation_results.SP[i] = u[i];
			}
			else if (simulation_type == h_cl) {
				simulation_results.P[i] = P_exc[nSamp + 1];
				simulation_results.OP[i] = OP_exc[nSamp + 1];
				simulation_results.SP[i] = simulation_results.SP[i - 1];
			}
			simulation_results.t[i] = i * Ts + t0;
			break;
		}


		// reinitialize intermediate vectors
		if (simulation_type == ol) {
			upsample(u[i], u[i + 1], nSamp, &P_us[0]);
			P_exc_old = P_exc[nSamp + 1];
			P_exc_old_old = P_exc[nSamp];
			if (i == 0) {
				P_exc[0] = u[0];
				P_exc[1] = u[0];
			}
			else {
				P_exc[0] = P_exc_old_old;
				P_exc[1] = P_exc_old;
			}
			for (size_t ct = 0; ct < nSamp; ct++) {
				P_exc[ct + 2] = P_us[ct];
			}
		}
		else if (simulation_type == cl) {
			real OP{ 0 };
			if (i == 0) {
				OP_exc[0] = 0;
				OP_exc[1] = 0;
			}
			else {
				OP = controller.pid(u[i], (simulation_results.x[i - 1] - x_min) / (x_max - x_min) * 100 + randn(gen_normal), i);
				if (OP > 100)
					OP = 100;
				else if (OP < 0)
					OP = 0;
				OP_exc[0] = simulation_results.OP[i - 1];
				OP_exc[1] = simulation_results.OP[i - 1];
			}
			for (size_t ct = 2; ct < nSamp + 2; ++ct) {
				OP_exc[ct] = OP;
				P[ct] = (dt * (p_max - p_min) / 100 * OP_exc[ct - 1] + tau_ip * P[ct - 1]) / (dt + tau_ip);
				P_exc[ct] = P[ct] + p_min;
			}
			P[0] = P[nSamp];
			P[1] = P[nSamp + 1];
			P_exc[0] = P_exc[nSamp];
			P_exc[1] = P_exc[nSamp + 1];
			simulation_results.OP[i] = OP;
			simulation_results.SP[i] = u[i];
		}
		else if (simulation_type == h_cl) {
			real SP{ 0.0 };
			if (i > 0) {
				SP = hydraulic_model(u[i], simulation_results.x[i - 1], i);
				simulation_results.x[i - 1] += randn(gen_normal) / 100 * (x_max - x_min);
			}
			else
				SP = 0;

			// Stem position controller
			real OP{ 0 };
			if (i == 0) {
				OP_exc[0] = 0;
				OP_exc[1] = 0;
			}
			else {
				OP = controller.pid(SP, (simulation_results.x[i - 1] - x_min) / (x_max - x_min) * 100, i);
				if (OP > 100)
					OP = 100;
				else if (OP < 0)
					OP = 0;
				OP_exc[0] = simulation_results.OP[i - 1];
				OP_exc[1] = simulation_results.OP[i - 1];
			}
			for (size_t ct = 2; ct < nSamp + 2; ++ct) {
				OP_exc[ct] = OP;
				P[ct] = (dt * (p_max - p_min) / 100 * OP_exc[ct - 1] + tau_ip * P[ct - 1]) / (dt + tau_ip);
				P_exc[ct] = P[ct] + p_min;
			}
			P[0] = P[nSamp];
			P[1] = P[nSamp + 1];
			P_exc[0] = P_exc[nSamp];
			P_exc[1] = P_exc[nSamp + 1];
			simulation_results.OP[i] = OP;
			simulation_results.SP[i] = SP;
		}


		x[0] = x[nSamp];
		x[1] = x[nSamp + 1];
		v[0] = v[nSamp];
		v[1] = v[nSamp + 1];
		a[0] = a[nSamp];
		a[1] = a[nSamp + 1];
		F_at[0] = F_at[nSamp];
		F_at[1] = F_at[nSamp + 1];
		F_res[0] = F_res[nSamp];
		F_res[1] = F_res[nSamp + 1];
		z_1[0] = z_1[nSamp];
		z_1[1] = z_1[nSamp + 1];

		simulation_results.x[i] = x[1];
		simulation_results.P[i] = P_exc[1];
		simulation_results.t[i] = i * Ts + t0;

		if (isnan(simulation_results.x[i])) {
			delete[] x;
			delete[] P_exc;
			delete[] P_us;
			delete[] v;
			delete[] a;
			delete[] F_at;
			delete[] F_res;
			delete[] z_1;
			delete[] OP_exc;
			delete[] P;
			clear_sim_data();
			return;
		}

		for (size_t j = 2; j < nSamp + 2; j++) {
			if (v[j - 1] > 0)
				sinal = 1;
			else if (v[j - 1] < 0)
				sinal = -1;
			else
				sinal = 0;

			if (v[j - 2] == 0.0 && stick_1 == 0)
				stick_1 = false;
			else if (v[j - 2] * v[j - 1] <= 0)
				stick_1 = true;

			real v_vs = v[j - 1] / v_s;
			real exp_v_s = std::exp(-v_vs * v_vs);
			s_v = (F_c + Fs_m_Fc * exp_v_s);

			// First element
			if (stick_1) {
				dot_z_1 = v[j - 1];
				z_1[j] = dt_u_2 * (v[j - 1] + v[j - 2]) + z_1[j - 1];
			}
			else {
				dot_z_1 = C_u_kappa1 * (sinal - F_at_1_ant / s_v);
				z_1[j] = dt_u_2 * (dot_z_1 + dot_z_1_ant) + z_1[j - 1];
			}
			F_at_1 = kappa_1 * z_1[j] + nu_1 * dot_z_1;
			if (std::abs(F_at_1) > s_v)
				stick_1 = false;

			F_at[j] = F_at_1 + F_v * v[j - 1];

			F_res[j] = S_a * P_exc[j - 1] - k * x[j - 1] - F_init - F_at[j];

			a[j] = inv_m * F_res[j];
			v[j] = dt_u_2 * (a[j] + a[j - 1]) + v[j - 1];
			x[j] = dt_u_2 * (v[j] + v[j - 1]) + x[j - 1];

			if (x[j] < x_min) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_min;
			}
			else if (x[j] > x_max) {
				F_res[j] = 0;
				a[j] = 0;
				v[j] = 0;
				x[j] = x_max;
			}

			dot_z_1_ant = dot_z_1;
			F_at_1_ant = F_at_1;
		}

	}
	delete[] x;
	delete[] P_exc;
	delete[] P_us;
	delete[] v;
	delete[] a;
	delete[] F_at;
	delete[] F_res;
	delete[] z_1;
	delete[] OP_exc;
	delete[] P;
}


template<typename real>
void ValveModel<real>::clear_sim_data() {
	simulation_results.t.clear();
	simulation_results.OP.clear();
	simulation_results.P.clear();
	simulation_results.x.clear();
	simulation_results.SP.clear();
	simulation_results.t.shrink_to_fit();
	simulation_results.OP.shrink_to_fit();
	simulation_results.P.shrink_to_fit();
	simulation_results.x.shrink_to_fit();
	simulation_results.SP.shrink_to_fit();
	sim_data_initialized = false;
}

template<typename real>
void ValveModel<real>::allocate_sim_data(size_t len_u) {
	if (!sim_data_initialized) {
		simulation_results.t.reserve(len_u);
		simulation_results.OP.reserve(len_u);
		simulation_results.P.reserve(len_u);
		simulation_results.x.reserve(len_u);
		simulation_results.SP.reserve(len_u);
		for (size_t i = 0; i < len_u; i++) {
			simulation_results.t.push_back(0.0);
			simulation_results.OP.push_back(0.0);
			simulation_results.P.push_back(0.0);
			simulation_results.x.push_back(0.0);
			simulation_results.SP.push_back(0.0);
		}
		sim_data_initialized = true;
	}
}


template<typename real>
std::vector<real> ValveModel<real>::OP2P_1order(std::vector<real>* OP) {

	std::vector<real> P;

	P = filter1order<real>(OP, this->get_tauip(), this->Ts);
	for (size_t i = 0; i < P.size(); ++i) {
		P[i] = P[i] * (p_max - p_min) / 100 + p_min;
		if (P[i] < p_min)
			P[i] = p_min;
		else if (P[i] > p_max)
			P[i] = p_max;
	}
	return P;
}


template<typename real>
std::vector<real> ValveModel<real>::filter2orderZP(const std::vector<real>* data, real wn, real xi) {

	real a0 = 1 + 2 * xi * wn * Ts + pow(Ts, 2) * pow(wn, 2);
	real b0 = (pow(Ts, 2) * pow(wn, 2)) / a0;
	real a1 = (2 + 2 * xi * wn * Ts) / a0;
	real a2 = -1 / a0;

	std::vector<real> filt_data((*data).size(), 0.0);
	filt_data[0] = ((*data)[0]);
	filt_data[1] = ((*data)[1]);

	for (size_t i = 2; i < data->size(); ++i)
		filt_data[i] = b0 * (*data)[i - 1] + a1 * filt_data[i - 1] + a2 * filt_data[i - 2];

	std::vector<real> inv_data(filt_data.size(), 0.0), filt2(filt_data.size(), 0.0);
	for (size_t i = 0; i < filt_data.size(); ++i)
		inv_data[i] = filt_data[filt_data.size() - 1 - i];

	filt2[0] = filt_data[0];
	filt2[1] = filt_data[1];

	for (size_t i = 2; i < data->size(); ++i)
		filt2[i] = b0 * inv_data[i - 1] + a1 * filt2[i - 1] + a2 * filt2[i - 2];

	std::vector<real> retdata(filt2.size(), 0.0);
	for (size_t i = 0; i < filt2.size(); ++i)
		retdata[i] = filt2[filt2.size() - 1 - i];

	return retdata;
}

template<typename real>
std::vector<real> ValveModel<real>::kalman_filter(const std::vector<real>* u, const std::vector<real>* y, const real Rv, const real Rw) {

	std::vector<real> P_norm(y->size(), 0.0);
	// Normalize the diaphragm pressure
	for (size_t i = 0; i < y->size(); ++i)
		P_norm[i] = (*y)[i] - p_min;

	real a = std::exp(-1 / tau_ip * Ts);
	real b = (p_max - p_min) / 100 * Ts / tau_ip;
	real c = 1;

	std::vector<real> M(u->size(), 0.0), P(u->size(), 0.0), Lc(u->size(), 0.0), hat_x(u->size(), 0.0), bar_x(u->size(), 0.0);
	M[0] = (real) 1.0;
	P[0] = (real) 1.0;
	Lc[0] = (real) 0.0;
	hat_x[0] = P_norm[0];
	bar_x[0] = P_norm[0];

	for (size_t i = 1; i < u->size(); ++i) {
		P[i] = M[i - 1] - M[i - 1] * c / (c * M[i - 1] * c + Rv) * c * M[i - 1];
		Lc[i] = P[i] * c / Rv;
		hat_x[i] = bar_x[i - 1] + Lc[i] * (P_norm[i - 1] - c * bar_x[i - 1]);
		bar_x[i] = a * hat_x[i] + b * (*u)[i - 1];
		M[i] = a * P[i] * a + Rw;
	}

	std::vector<real> retdata(y->size(), 0.0);
	for (size_t i = 0; i < y->size(); ++i)
		retdata[i] = hat_x[i] + p_min;

	return retdata;
}


template<typename real>
real ValveModel<real>::hydraulic_model(real SP, real x, size_t ct) {
	int seed = std::chrono::system_clock::now().time_since_epoch().count() + 1213155;
	std::default_random_engine gen_normal(seed);
	std::normal_distribution<real> randn(0.0, 50 / pow(10, 25 / 10));

	if (ct == 0) {
		Q_int[0] = 100;
		Q[0] = 100;
		return 0;
	}
	else {
		real tau_h = 3.375;
		Q_int[ct] = -100442.87 * x * x - 630.57 * x + 101.54;
		Q[ct] = (Ts * Q_int[ct] + tau_h * Q[ct - 1]) / (tau_h + Ts);
		return controller_hydraulic.pid(SP, Q[ct] + randn(gen_normal), ct);
	}
}