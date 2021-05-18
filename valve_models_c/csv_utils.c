#include "stdio.h"
#include "stdlib.h"
#include "csv_utils.h"

csv_data get_data(const char* filename, int n){
    

    float *input = calloc(n, sizeof(float));
    int count = 0;

    float t, OP, x_kano, x_karnopp, x_lugre, x_gms;

    FILE *fp;
    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error reading file\n");
    }

    char line[2048];
    while (!feof(fp)){
        fgets(line, 2048, fp);
        sscanf(line, "%e,%e,%e,%e,%e,%e,%e\n", &t, &OP, &input[count], &x_kano, &x_karnopp, &x_lugre, &x_gms);
        count++;
    }

    csv_data obtained_data;
    obtained_data.data = calloc(count, sizeof(double));
    obtained_data.n = count;
    for (int i=0; i<count; i++){
        obtained_data.data[i] = (double) input[i];
    }
    
    fclose(fp);
    free(input);

    return obtained_data;
}