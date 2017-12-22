//util test
// file: util_test.cc
//

// isip include files
//
#include "asc_utils.h"
#include <limits.h>

// spectrogram: automated sound comparison utility
//
// This is a driver program that demonstrates the functionality
// of the sound comparison module.
//


int main(int argc, const char** argv) {

 
double **M=NULL;
double **N=NULL;
double **result=NULL;

double *p;
double *q;
long error_code=0;

p=(double *)malloc(10*sizeof(double));
for (int i=0;i<10;i++)
	p[i]=i+1;

printf ("p1: %p %p\n",&p,p);
printf ("q1: %p %p\n",&q,q);	
Asc_utils::array_sqrt(p,10,&q);


printf ("p2: %p %p\n",&p,p);
printf ("q2: %p %p\n",&q,q);

Asc_utils::array_sqrt(p,10,&p);	
	
printf ("p3: %p %p\n",&p,p);
printf ("q3: %p %p\n",&q,q);


free(p);

printf("*");

free(q);

//allocate memory for D
/* M = (double **)malloc( 5* sizeof(double *));
for (int i = 0; i < 5; i++) 
M[i] = (double *)malloc(5 * sizeof(double));
 
for (int i=0;i<5;i++)
	for (int j=0;j<5;j++)
	{
		if (i==j)
		  M[i][j]=i+1;
		else  
		  M[i][j]=0;
	}
M[0] [1]=1;

N = (double **)malloc( 5* sizeof(double *));
for (int i = 0; i < 5; i++) 
N[i] = (double *)malloc(3 * sizeof(double));
 
for (int i=0;i<5;i++)
	for (int j=0;j<3;j++)
	{
		if (i==j)
		  N[i][j]=i-1;
		else  
		  N[i][j]=0;
	}
N[0] [1]=1;
if ((error_code=Asc_utils::matrix_multiply(M,5,5,N,5,3,&result)) !=0)
{
fprintf(stdout,"Error : %ld\n",error_code);
exit(-1);
}   */
 
/* if ((error_code=Asc_utils::DP_path(&p,&q,M,5,5)) !=0)
{
fprintf(stdout,"Error : %ld\n",error_code);
exit(-1);
}

 for (int i=0;i<5;i++)
	printf("\np[%d]=%d\t",i,p[i]);

for (int i=0;i<5;i++)
	printf("\nq[%d]=%d\t",i,q[i]);  */
	
//free(p);
//free(q);	
/* double test_v;

printf("%f\n",test_v);

N[0][0]=0*1.0/0;
printf("N before %p %p\n",N,*N);
Asc_utils::matrix_setNaN(N,5,3,100,&N);
printf("\nN after %p %p\n",N,*N);

for (int i=0;i<5;i++)
	{
	for (int j=0;j<5;j++)
	 printf("%f\t",M[i][j]);
	printf("\n");
	}
	printf("\n");
for (int i=0;i<5;i++)
	{
	for (int j=0;j<3;j++)
	 printf("%f\t",N[i][j]);
	printf("\n");
	}	
printf("\n"); */	


	
/* for (int i=0;i<5;i++)
	{
	for (int j=0;j<3;j++)
	 printf("%f\t",result[i][j]);
	printf("\n");
	} */
	 


/* for (int k = 0; k < 5; k++) 
	free(M[k]);
  free(M);
  
  for (int k = 0; k < 5; k++) 
	free(N[k]);
  free(N);
  
  for (int k = 0; k < 5; k++) 
	free(result[k]);
  free(result);  */
  
 }
