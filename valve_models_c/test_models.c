#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "time.h"


#include <unistd.h>


#include "valve_models.h"
#include "csv_utils.h"


#define PARAM_VALVULA {1.6, 210490, 0.0445, 2550, 0, 0.029, 180000, 41368, 3.571428571428571}
#define PARAM_ATRITO {700, 780, 125000, 5.0e-04}

#define TS 1e-3

int main(int argc, char* argv[]) {

    char cwd[1024];
    getcwd(cwd, sizeof(cwd));
    chdir(strcat(cwd, "/test_data"));

    csv_data dados;
	const char* file = "dados.csv";
	int n = 3600/TS;
	dados = get_data(file, n);

    double param_valvula[] = PARAM_VALVULA;
    double param_atrito[] = PARAM_ATRITO;
    double pos0[] = {0, 0};
    double dt = 1e-5;

    clock_t t; 
    t = clock(); 
    results *simdata = sim_karnopp(dados.data, dados.n, (double) TS, param_valvula, param_atrito, pos0, dt);
    t = clock() - t; 
    double time_taken = ((double)t)/CLOCKS_PER_SEC;
    printf("exec time: %e\n", time_taken);

    FILE *csv_file;
    csv_file = fopen("simulation.csv", "w+");
    for (int i = 0; i < dados.n; i++){
		fprintf(csv_file, "%e,%e,%e,%e,%e,%e\n", simdata->t[i], simdata->P[i], simdata->x[i], simdata->v[i], simdata->a[i], simdata->F_at[i]);
	}
    fclose(csv_file);
	return 0;
}