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

#include "functions.h"

void prepararObjFunc(int* FUNCTION, double* lb, double* ub)/*{{{*/
{
	switch (*FUNCTION) {
		case 0: //Rastrigin
			*lb = -5.12;
			*ub = 5.12;
			break;
		case 1: //Schaffer
			*lb = -100.00;
			*ub = 100.00;
			break;
		case 2: //Griewank
			*lb = -600.00;
			*ub = 600.00;
			break;
		case 3: //Ackley
			/*
			   -	Dimension: n arbitrary
			   -       Domain:   -32 <= | x_i | <= 32.0
			   -       Minimum 0 at x_i = 0.0
			   */
			*lb = -32.00;
			*ub = 32.00;
			break;
		case 4: //Rosenbrock
			*lb = -30.00;
			*ub = 30.00;
			break;
		case 5: //Sphere
			*lb = -100.00;
			*ub = 100.00;
			break;
		case 6: //MPE
			*lb = 0;
			*ub = 5;
			break;
		case 7: //SCHAFFER_F6
			*lb = -100.00;
			*ub = 100.00;
			break;
		case 8: //G_SCHWEFELS
			*lb = -500.00;
			*ub = 500.00;
			break;
		case 12: // Levy Function
			*lb = -10.00;
			*ub = 10.00;
			break;
		case 13: // Zakharov
			*lb = -5.00;
			*ub = 10.00;
			break;
		case 14: // EggHolder
			*lb = -512.00;
			*ub = 512.00;
			break;
		case 15: // Holzman
			*lb = -10.00;
			*ub = 10.00;
			break;
		case 16: // Michalewitz
			*lb = 0.00;
			*ub = PI;
			break;
		case 18: // Powell
			*lb = -4.00;
			*ub = 5.00;
			break;
		case 19: // Rana
			*lb = -512.00;
			*ub = 512.00;
			break;
		case 20: // Shubert
			*lb = -10.00;
			*ub = 10.00;
			break;
		case 77: // Schwefels 222
			*lb = -10.00;
			*ub = 10.00;
			break;
		default:
			printf("Info: Invalid function.\n") ;
			exit(0);
	}
}/*}}}*/

double objfunc(double sol[],const int* FUNCTION, const int* DIM, int *cont)/*{{{*/
{
	*cont+=1;

	switch (*FUNCTION) {
		case 0: //Rastrigin

			return rastrigin(sol,*DIM);

		case 1: //Schaffer

			return schaffer(sol,*DIM);

		case 2: //Griewank

			return griewank(sol,*DIM);

		case 3: //Ackley

			return ackley(sol,*DIM);

		case 4: //Rosenbrock

			return rosenbrock(sol,*DIM);

		case 5: //Sphere

			return sphere(sol,*DIM);

		case 6: //MPE
			return mpe(sol,*DIM);

		case 7: //SCHAFFER_F6

			return schaffer_f6(sol,*DIM);

			break;
		case 8: //G_SCHWEFELS

			return g_schwefels(sol,*DIM);

			break;

		case 12: // Levy Function
			
			return levy(sol,*DIM);
			break;
		case 13: // Zakharov
			return zakharov(sol,*DIM);
			break;
		case 14: // EggHolder

			return egg_holder(sol,*DIM);
			break;
		case 15: // Holzman
			
			return holzman(sol,*DIM);
			break;

		case 16: // Michalewitz
			
			return michalewitz(sol,*DIM);
			break;

		case 18: // Powell
			return powell(sol,*DIM);
			break;

		case 19: // Rana

			return rana(sol,*DIM);
			break;
		case 20: // Shubert

			return shubert(sol,*DIM);
			break;
		case 77: //  Schwefel's function 2.22

			return schwefels222(sol,*DIM);
			break;
	
		default:
			printf("Info: Invalid function.\n") ;
			exit(0);
	}
}/*}}}*/

double rastrigin( double SOL[],  int DIM){/*{{{*/
	//SOL[] -> the population set
	//DIM the dimension of SOL[]
	int j;
	double top = 0.00;

	for(j=0;j<DIM;j++)
	{
		top=top+(pow(SOL[j],(double)2)-10*cos(2*M_PI*SOL[j])+10);
	}
	return top;
}/*}}}*/

double schaffer( double SOL[],  int DIM){/*{{{*/
	//SOL[] -> the population set
	//DIM the dimension of SOL[]

	int j;
	double top = 0.00 , top1 = 0.00;

	for(j=0;j<DIM;j++)
	{
		top=top+(pow(SOL[j],(double)2));
	}

	top = pow(top,(double)0.25);

	for(j=0;j<DIM;j++)
	{
		top1=top1+(pow(SOL[j],(double)2));
	}

	top1=pow(top1,(double)0.1);
	top1 = pow(sin(50*top1),(double)2) +1.0;

	return top*top1;
}/*}}}*/

double griewank( double SOL[],  int DIM){/*{{{*/
	//SOL[] -> the population set
	//DIM the dimension of SOL[]

	int j;
	double top = 0.00 , top1 = 0.00, top2 = 0.00;

	for(j=0;j<DIM;j++)
	{
		top1=top1+pow((SOL[j]),(double)2);
		top2=top2*cos((((SOL[j])/sqrt((double)(j+1)))*M_PI)/180);
	}
	top=(1/(double)4000)*top1-top2+1;

	return top;
}/*}}}*/

double ackley( double SOL[],  int DIM){/*{{{*/
	//SOL[] -> the population set
	//DIM the dimension of SOL[]

	int i;	
	double aux = 0.0;
	double aux1 = 0.0;

	for (i = 0; i < DIM; i++)
	{
		aux += SOL[i]*SOL[i];
	}
	for (i = 0; i < DIM; i++)
	{
		aux1 += cos(2.0*M_PI*SOL[i]);
	}

	return (-20.0*(exp(-0.2*sqrt(1.0/(float)DIM*aux)))-exp(1.0/(float)DIM*aux1)+20.0+exp(1));
}/*}}}*/

double rosenbrock( double SOL[],  int DIM){/*{{{*/
	//SOL[] -> the population set
	//DIM the dimension of SOL[]

	int i;
	double top = 0.00;

	for (i = 0; i < DIM-1; i++)
	{
		top=top+100.*pow((SOL[i+1] - pow(SOL[i],2.)),2) + pow((1. - SOL[i]),2);
	}

	return top;
}/*}}}*/

double sphere( double SOL[],  int DIM){/*{{{*/
	//SOL[] -> the population set
	//DIM the dimension of SOL[]

	int j;
	double top = 0.00;

	for(j=0;j<DIM;j++)
	{
		top=top+SOL[j]*SOL[j];
	}

	return top;
}/*}}}*/

double shubert(double SOL[],  int DIM){/*{{{*/
	//SOL[] -> the population set
	//DIM the dimension of SOL[]
	//Shubert
	/*
	   -   Domain  |x| <= 10.0
	   -   Number of local minimum = 400
	   -   Global minimum fmin = -24.062499 at the ff. points
	   -    (-6.774576, -6.774576), ..., (5.791794, 5.791794)
	   */
	double sum = 0.0;
	int i;
	for (i = 0; i < DIM; i++) {
		sum += -sin(2.0*SOL[i]+1.0)
			-2.0*sin(3.0*SOL[i]+2.0)
			-3.0*sin(4.0*SOL[i]+3.0)
			-4.0*sin(5.0*SOL[i]+4.0)
			-5.0*sin(6.0*SOL[i]+5.0);
	}

	return sum/(DIM/2.0);
}/*}}}*/

char *getFunctionName(int FUNCTION){/*{{{*/

	switch (FUNCTION) {
		case 0: //Rastrigin

			return "Rastrigin";

		case 1: //Schaffer

			return "Schaffer F7";

		case 2: //Griewank

			return "Griewank";

		case 3: //Ackley

			return "Ackley";

		case 4: //Rosenbrock

			return "Rosenbrock";

		case 5: //Sphere

			return "Sphere";

		case 6: //MPE
			return "Molecular Potential Energy";
		case 7: //SCHAFFER_F6

			return "Schaffer F6";

			break;
		case 8: //G_SCHWEFELS

			return "Generalized Schwefels";

			break;
		case 20: // Shubert

			return "Shubert";
			break;
		case 16: // Michalewitz

			return "Michalewitz";
			break;
		case 18: // Powell
			return "Powell";
			break;

		case 12: // Levy
			return "Levy";
			break;

		case 13: // Zakharov
			return "Zakharov";
			break;

		case 14: // Egg holder
			return "Egg holder";
			break;

		case 15: // Generalized Holzman
			return "Generalized Holzman";
			break;

		case 19: // Rana
			return "Rana";
			break;
		case 77: //  Schwefel's function 2.22
			return " Schwefel's function 2.22";
			break;


		default:
			printf("Info: Invalid function.\n") ;
			exit(0);
	}
}/*}}}*/

double schaffer_f6( double sol[],  int DIM){/*{{{*/
	double top=0;
	double top1=0;
	double top2=0;
	int j;
	for(j=0;j<DIM-1;j++)
	{
		top1 = pow( sin( sqrt( pow(sol[j],(double)2) + pow(sol[j+1],(double)2) ) ), (double)2) - 0.5;
		top2 = pow( 1.0 + 0.001* ( pow(sol[j],(double)2) + pow(sol[j+1],(double)2) ), (double)2);
		top = top + (0.5 + top1/top2);
	}
	return top;
}/*}}}*/

double g_schwefels( double sol[],  int DIM){/*{{{*/
	//known_optimal = -418.982887272433 at sol(i)=420.9687
	int i;
	double aux = 0.0;
	for (i=0;i<DIM;i++)
	{
		aux += (-1 * sol[i]) * sin(sqrt(fabs(sol[i])));
	}
	return(aux);
}/*}}}*/

double mpe( double sol[], int DIM){/*{{{*/
	//[0,5] fmin=-0.411183034 * n
	int i;
	double aux1=0.0,aux2=0.0;
	for (i=0;i<DIM;i++)
	{
		aux1 += 1 + cos(3*sol[i]);
		aux2 += (pow(-1,i) / (sqrt(fabs(10.60099896-4.141720682 * cos(sol[i])))));
	}
	return aux1 + aux2;
}/*}}}*/

double michalewitz( double sol[], int DIM){/*{{{*/
	int i;
	double aux;
	//Michalewitz
	/*
Dimension: n
Domain: 0< | x[i] | <= PI
Global minimum: x[] = -0.966*n
*/
	aux=0;
	for (i=0;i<DIM;i++) {
		aux = aux + sin(sol[i]) * pow(sin((i+1)*sol[i]*sol[i]/(float)PI), 2.0 * 10.0);
	}
	return(-1*aux);
}/*}}}*/

double powell( double sol[], int DIM){/*{{{*/
	int j;
	double aux; //Powell
	/*
Dimension: n > 4
Domain: -4<= x[i] <= 5
Global minimum: at (3, -1, 0, 1, ..., 3, -1, 0, 1) with fmin = 0
*/
	aux = 0.0;
	for (j = 1; j <= (int)DIM/4; j++) {
		aux +=  pow(sol[4*j-4] + 10 * sol[4*j-3],2.0)
			+ 5 * pow(sol[4*j-2] - sol[4*j-1],2.0)
			+ pow(sol[4*j-3] - 2 * sol[4*j-2], 4.0)
			+ 10 * pow(sol[4*j - 4] - sol[4*j-1], 4.0);
	}
	return aux;
}/*}}}*/

double levy( double sol[], int DIM)/*{{{*/
{
	//Levy Function
	double aux,*y;
	int i;
	//x[i] = 1 f(x[i])=0
	y = (double*) malloc (DIM * sizeof(double));
	for (i = 0; i< DIM; i++)
		y[i] = 1+(sol[i]-1)/4.0;
	aux = pow(sin(PI*y[0]),2.0);
	for (i = 0; i<DIM-1;i++)
		aux = aux + pow(y[i]-1,2.0)*(1+10*pow(sin(PI*y[i]+1),2.0));

	aux = aux+pow(y[DIM-1]-1,2.0)*(1+pow(sin(2*PI*y[DIM-1]),2.0) );

	free (y);
	return ( aux );
}/*}}}*/

double zakharov( double sol[], int DIM)/*{{{*/
{
	// Zakharov function  //x[i] = 0 f(x[i])=0
	//
	double aux,aux1;
	int j;
	aux = aux1 = 0.0;
	for (j = 0; j< DIM; j++)
	{
		aux = aux + pow(sol[j],2.0);
		aux1 = aux1+0.5*j*sol[j];
	}


	return ( aux+pow(aux1,2.0)+pow(aux1,4.0) );
}/*}}}*/

double egg_holder( double sol[], int DIM)/*{{{*/
{
	//Egg holder
	/*
	   - Dimension: n
	   - Domain:  -512 < | x_i | < 512
	   - Minimum for n=2 fmin(512, 404.2319) = -959.641
	   */
	double aux;
	int i;
	aux = 0.0;
	for (i = 0; i < DIM-1; i++) 
	{
		aux += -(sol[i+1] + 47.0) * sin(sqrt(fabs(sol[i+1] + sol[i] * 0.5 + 47.0))) + sin(sqrt(fabs(sol[i] - (sol[i+1] + 47.0)))) * (-sol[i]);
	}
	return (aux);
}/*}}}*/

double rana( double sol[], int DIM)/*{{{*/
{
	//Rana
	/*
Dimension: n
Domain: -520<= x[i] <= 520
Global minimum: ???
*/
	double sum,t1,t2;
	int i;
	sum = 0.0;
	for (i = 0; i < DIM-1; i++) {
		/* if(sol[i] < -520.00f || sol[i] > 520.00f) */
		/* 	printf("rana %f\n",sol[i]); */

		t1 = sqrt(fabsf(sol[i+1] + sol[i] + 1.0));
		t2 = sqrt(fabsf(sol[i+1] - sol[i] + 1.0));
		sum += (sol[i+1] + 1.0) * cos(t2) * sin(t1) + cos(t1) * sin(t2) * sol[i];
	}
	return sum/(double)(DIM-1);
}/*}}}*/

double holzman( double sol[], int DIM){/*{{{*/
	//Generalized Holzman
	/*
	   Dimension: n
	   Domain: | x[i] | <= 10
	   Global minimum: 0 at x[i] = 0
	   */
	int i;
	double aux = 0.0;
	for (i = 0; i < DIM; i++) 
	{
	aux += i * pow(sol[i] , 4);
	}
	return aux;
}/*}}}*/

double schwefels222( double sol[], int DIM){/*{{{*/
	// Schwefel's function 2.22
	/*
	   - Domain:  | x_i | <= 10.0
	   Global minimum is 0.0 at x_i = 0.00
	   */

	int i;
	double aux = 0.0,aux1=0.0;
	for (i=0;i<DIM;i++)
	{
		aux += fabs(sol[i]);
		aux1 *= fabs(sol[i]);
	}

	return (aux+aux1);
}/*}}}*/

