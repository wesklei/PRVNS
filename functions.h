/* 
        Description: The ANSI C code of the PRVNS approach
        Programmer:  Wesklei Migliorini
        E-mail:      wesklei.m@gmail.com
        Date:	     04/11/2014
        Lisence:     Free
        Note:        The system was developed using Linux.
        To compile:  Type: make
        To run: ./algorithm input.in
 */
#ifndef functions_h
#define functions_h

#include "data.h"
#define PI 3.14159265

void prepararObjFunc(int* FUNCTION, double* lb, double* ub);
double objfunc(double sol[], const int* FUNCTION, const int* DIM, int *cont);
char *getFunctionName(int FUNCTION);
double rastrigin( double SOL[],  int DIM);
double schaffer( double SOL[],  int DIM);
double griewank( double SOL[],  int DIM);
double ackley( double SOL[],  int DIM);
double rosenbrock( double SOL[],  int DIM);
double sphere( double SOL[],  int DIM);
double schaffer_f6( double sol[],  int DIM);
double g_schwefels( double sol[],  int DIM);
double mpe( double sol[],  int DIM);
double shubert(double SOL[],  int DIM);
double powell( double sol[], int DIM);
double michalewitz( double sol[], int DIM);
double levy( double sol[], int DIM);
double zakharov( double sol[], int DIM);
double egg_holder( double sol[], int DIM);
double generalized_holzman( double sol[], int DIM);
double rana( double sol[], int DIM);
double holzman( double sol[], int DIM);
double schwefels222( double sol[], int DIM);
double shifted_sphere( double sol[], int DIM);
double shifted_schwefel_221( double sol[], int DIM);
double shifted_rosenbrock( double sol[], int DIM);
double shifted_rastrigin( double sol[], int DIM);
double shifted_griewank( double sol[], int DIM);
double shifted_ackley( double sol[], int DIM);
double shifted_schwefel_222(double sol[], int dim); 
double shifted_schwefel_12(double sol[], int dim);
double f_10(double x, double y);// auxiliar
double shifted_extended_f10(double sol[], int dim);
double shifted_bohachevsky(double sol[], int dim);
double shifted_schaffer(double sol[], int dim);


#endif
