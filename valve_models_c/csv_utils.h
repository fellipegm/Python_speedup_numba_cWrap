#pragma once
#ifndef _CSV_UTILS_H
#define _CSV_UTILS_H
typedef struct csv_data
{
    double* data;
    int n;
} csv_data;


csv_data get_data(const char* filename, int n);

#endif //_CSV_UTILS_H