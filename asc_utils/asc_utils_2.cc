// file: asc_utils_2.cc
//
// This file contains implementations of the private methods.
//


#include <string.h>

// isip include files
//
#include "asc_utils.h"


//------------------------------------------------------------------------
//
// public methods
//
//-----------------------------------------------------------------------

// method: DP_path
// 
// arguments:
//  int *p,*q: output state sequences  
// Originally from Matlab code : Copyright (c) 2003 Dan Ellis <dpwe@ee.columbia.edu>
// return: a value of Asc_utils::ERR if there is an error
// return mn if no error
//Use dynamic programming to find min cost path through matrix M (m by n)
 long Asc_utils::DP_path(long **p,long **q,double **M,int m,int n)
{

long mn=max(m,n);
//if (*p==NULL)
//*p = (long *)malloc(m*n* sizeof(long));
//if(*q==NULL)
//*q = (long *)malloc(m*n* sizeof(long));

double tmp[3]; // temp array

long tb;
double dmax;

if (*p == NULL || *q == NULL)
{
	return Asc_utils::ERR_MEM;
}	
double **D;
//allocate memory for D
D = (double **)malloc((m+1)* sizeof(double *));
for (int i = 0; i < (m+1); i++) 
D[i] = (double *)malloc((n+1) * sizeof(double));
//setting D to zero
//for(int i=0;i<(m+1);i++)
	//for(int j=0;j<(n+1);j++)
	// D[i][j]=0;
for(int j=0;j<(n+1);j++)
	 D[0][j]=DBL_MAX;	 
for(int i=0;i<(m+1);i++)
	 D[i][0]=DBL_MAX;	
D[0][0]=0;
 for(int i=1;i<(m+1);i++)
	for(int j=1;j<(n+1);j++)
	  D[i][j]=M[i-1][j-1];
	  
//trace back
int **phi;
//allocate memory for D
phi = (int **)malloc((m)* sizeof(int *));
for (int i = 0; i < (m); i++) 
phi[i] = (int *)malloc((n) * sizeof(int));	  


if (phi == NULL || D == NULL)
{
	return Asc_utils::ERR_MEM;
}


//setting phi to zero
for(int i=0;i<(m);i++)
	for(int j=0;j<(n);j++)
		phi[i][j]=-1;
for(int i=0;i<(m);i++)
{
	for(int j=0;j<(n);j++)	
    {	
		tmp[0]=D[i][j];
		tmp[1]=D[i][j+1];
		tmp[2]=D[i+1][j];
		tb=array_min(tmp,3);
		
        
		dmax=tmp[tb];
	    
		
		
		D[i+1][j+1]=D[i+1][j+1]+dmax;
		phi[i][j]=tb;
	}
}



long i=m-1;
long j=n-1;
long count=0;

(*p)[count]=i;
(*q)[count]=j;
count++;
while( i>0 && j>0)
{

	tb=phi[i][j];
	
	if (tb == 0)
	{
		i--; 
		j--;
	}
	else if (tb == 1)
	{
		i--;
	}
	else if (tb == 2)
	{ 
		j--;
	}
	else
	{
		return Asc_utils::ERR;
	}
   

	
(*p)[count]=i;
(*q)[count]=j;
	count++;
}

//free 2-d arrays	
  for (int k = 0; k < m+1; k++) 
	free(D[k]);
  free(D);
  for (int k = 0; k < m; k++) 
	free(phi[k]);
  free(phi);

return count;  
}


// method: array_min
// 
// arguments: input array
// return:  index of min
//array_min
long Asc_utils::array_min(double *input,long len)
{
double min_t=input[0];
long index=0;
for(long i=0;i<len;i++)
{
	if (input[i] < min_t)
	{
		min_t=input[i];
		
		index=i;
	}
}

return index;
}

// method: array_sqrt
// 
// arguments: input array
// return:  result
//array_sqrt
long Asc_utils::array_sqrt(double *input,long len,double ** result)
{
	// (*result)=(double*)malloc(len*sizeof(double));
	for (int i=0;i<len;i++)
		(*result)[i]=sqrt(input[i]); 

	return Asc_utils::NO_ERR;  
}

long Asc_utils::array_pow2(double *input,long len,double ** result)
{
	// (*result)=(double*)malloc(len*sizeof(double));
	for (int i=0;i<len;i++)
		(*result)[i]=(input[i]*input[i]); 

	return Asc_utils::NO_ERR;  
}

// method: matrix_multiply
// 
// arguments:
//  input1 input2 and their sizes
//   output result
// Originally from Matlab code : Copyright (c) 2003 Dan Ellis <dpwe@ee.columbia.edu>
// return: a value of Asc_utils::ERR if there is an error
//
long Asc_utils::matrix_multiply(double **input1,int rows1,int cols1,double **input2,int rows2,int cols2,double ***result)
{
	double tmp;
	//if (result==NULL){
		*result = (double **)malloc(rows1 * sizeof(double *));
		for(int i = 0; i < rows1; i++)
		{
			(*result)[i] = (double *)malloc(cols2 * sizeof(double));
		}
		if (*result == NULL)
		{
			return  Asc_utils::ERR_MEM;
		}
	//}
	if( cols1 != rows2 )
	{

		fprintf(stdout,"Inputs sizes are not consistent.\n");
		return Asc_utils::ERR;
	}
    else
	{
		 for (int i=0;i<rows1;i++)
        {

            for(int j=0;j<cols2;j++)
            {
				tmp=0;
				for(int k=0;k<rows2;k++)
                {
                  
                    tmp += input1[i][k]*input2[k][j];

                }
			    (*result)[i][j]=tmp;
			
			}
		}
	
	}
	
return Asc_utils::NO_ERR; 	
}   


// method: matrix_sum
// 
// arguments:
//  input1, row1,col1,dim
//   output result
// 
// return: a value of Asc_utils::ERR if there is an error
// sum over one of the 2 dimesions of input matrix  dim=0 --> rows dim=1 --->columns
long Asc_utils::matrix_sum(double **input,int rows,int cols,int dim,double **result)
{
	double tmp;
	if (dim == 0)
	{
		//*result=(double *)malloc(cols* sizeof(double));
		for (int i=0;i<cols;i++)
		{
			tmp=0;
			for (int j=0;j<rows;j++)
				tmp+=input[j][i];
			(*result)[i]=tmp;	
		}		
	}
	else if (dim==1)
	{
		//*result=(double *)malloc(rows* sizeof(double));
		for (int i=0;i<rows;i++)
		{
			tmp=0;
			for (int j=0;j<cols;j++)
				tmp+=input[i][j];
			(*result)[i]=tmp;	
		}
	}
	else
	{
		return Asc_utils::ERR;
	}

return Asc_utils::NO_ERR;
}

// method: matrix_sum
// 
// arguments:
//  input1, rows,cols
//   output result
// 
// return: a value of Asc_utils::ERR if there is an error
// calculate the pow two of matrix elements
long Asc_utils::matrix_epow2(double **input,int rows,int cols,double ***result)
{
	
	//*result = (double **)malloc(rows * sizeof(double *));
	/* for(int i = 0; i < rows; i++)
	{
		(*result)[i] = (double *)malloc(cols * sizeof(double));
	} */
	for (int i=0;i<rows;i++)
		for(int j=0;j<cols;j++)
			(*result)[i][j]=input[i][j]*input[i][j];
			
	return Asc_utils::NO_ERR;
}

// method: matrix_setNaN
// 
// arguments:
//  input1, rows,cols ,set_to
//   output result
// 
// return: a value of Asc_utils::ERR if there is an error
// set NaN mambers to another value
long Asc_utils::matrix_setNaN(double **input,int rows,int cols,double set_to,double *** result)
{

  
	//*result = (double **)malloc(rows * sizeof(double *));
	for(int i = 0; i < rows; i++)
	{
		(*result)[i] = (double *)malloc(cols * sizeof(double));
	}
	
	
	
	
	for (int i=0;i<rows;i++)
		for(int j=0;j<cols;j++)
			if(isnan(input[i][j]) ==true)
				(*result)[i][j]=set_to;
			else
				(*result)[i][j]=input[i][j];
				
	return Asc_utils::NO_ERR;	
}


// method: matrix_copy
// 
// arguments:
//  
// 
// return: a value of Asc_utils::ERR if there is an error
// copy input to result
long Asc_utils::matrix_copy(double **input,int rows,int cols,double ***result)
{
	for(int i=0;i<rows;i++)
		for(int j=0;j<cols;j++)
			(*result)[i][j]=input[i][j];
    return Asc_utils::NO_ERR;
}


// method: array_copy
// 
// arguments:
//  
// 
// return: a value of Asc_utils::ERR if there is an error
// copy input to result
long Asc_utils::array_copy(double *input,int len,double **result)
{
	for(int i=0;i<len;i++)
			(*result)[i]=input[i];
    return Asc_utils::NO_ERR;
}

// method: matrix_sim
// 
// arguments: input1, input2 
//  
// Originally from matlab code : Copyright (c) 2003 Dan Ellis <dpwe@ee.columbia.edu>

// return: a value of Asc_utils::ERR if there is an error
// Calculate sim matrix  between input1 and input2

long Asc_utils::matrix_sim(double **input1,int rows1,int cols1,double **input2,int rows2,int cols2,double ***M)
{
        long error_code=0;
		double *EA;
		double *EB;
		double tmp;
		
		EA=(double*)malloc(cols1*sizeof(double));
		EB=(double*)malloc(cols2*sizeof(double));
		
		if (rows1 != rows2)
		{	
			return Asc_utils::ERR;
		}
	    /* *M = (double **)malloc(cols1 * sizeof(double *));
		for(int i = 0; i < cols1; i++)
		{
			(*M)[i] = (double *)malloc(cols2 * sizeof(double));
		} */
		
		double **tmp1;
		double *tmp2;
		tmp1 = (double **)malloc(rows1 * sizeof(double *));
	    for(int i = 0; i < rows1; i++)
	      tmp1[i] = (double *)malloc(cols1 * sizeof(double));
		tmp2=(double*)malloc(cols1*sizeof(double));
		
	   if ((error_code=matrix_epow2(input1,rows1,cols1,&tmp1)) !=Asc_utils::NO_ERR)
			return error_code;
	   
	   if((error_code=matrix_sum(tmp1,rows1,cols1,0,&tmp2)) !=Asc_utils::NO_ERR)
			return error_code;
	   
	   
	   if((error_code=array_sqrt(tmp2,cols1,&EA)) !=Asc_utils::NO_ERR)
			return error_code;
		
		
	   //free tmp vars
	   for (int i = 0; i < rows1; i++) 
			free(tmp1[i]);
	   free(tmp1);
	   free(tmp2);
		
        double **tmp3;
		double *tmp4;
		tmp3 = (double **)malloc(rows2 * sizeof(double *));
	    for(int i = 0; i < rows2; i++)
	      tmp3[i] = (double *)malloc(cols2 * sizeof(double));
		tmp4=(double*)malloc(cols2*sizeof(double));
	   if ((error_code=matrix_epow2(input2,rows2,cols2,&tmp3)) !=Asc_utils::NO_ERR)
			return error_code;
	  
	   if((error_code=matrix_sum(tmp3,rows2,cols2,0,&tmp4)) !=Asc_utils::NO_ERR)
			return error_code;
	   if((error_code=array_sqrt(tmp4,cols2,&EB)) !=Asc_utils::NO_ERR)
			return error_code;
			
	   //free tmp vars
	   for (int i = 0; i < rows2; i++) 
			free(tmp3[i]);
	   free(tmp3);
	   
	   free(tmp4);		
	 for (int i=0;i<cols1;i++)
		for(int j=0;j<cols2;j++)
		{
		   tmp=0;
		   for(int k=0;k<rows1;k++)
				tmp+=input1[k][i]*input2[k][j];
			if ((EA[i]*EB[j]) !=0) 	
				(*M)[i][j]=tmp/(EA[i]*EB[j]);
			else
				(*M)[i][j]=0; // zero out Nans
		}
		
	 free(EA);
     free(EB);	 
	 return Asc_utils::NO_ERR;	
}

// method: linear_corr
// 
// arguments: input1, input2 ,len
//  
// 

// return: a value of Asc_utils::ERR if there is an error
//linear correlation
long Asc_utils::linear_corr(double *x,double* y,int len,double *r)
{
	double mx=0;
	double my=0;
	double *tmp1=NULL,*tmp2=NULL,*tmp3=NULL;
	double tmp4,tmp5,tmp6;
	long count=0;
	
	tmp1=(double *)malloc(len* sizeof(double));
	tmp2=(double *)malloc(len* sizeof(double));
	tmp3=(double *)malloc(len* sizeof(double));
	
	for (int i=0;i<len;i++)
	{	
		if(isnan(x[i])==0 && isnan(y[i])==0)
		{
			mx+=x[i];
			my+=y[i];
			count++;
		}
    }
	mx=mx/(double)count;
	my=my/(double)count;
	count=0;
    for (int i=0;i<len;i++)
	{
		if(isnan(x[i])==0 && isnan(y[i])==0)
		{
			tmp1[count]=(x[i]-mx)*(y[i]-my);
			tmp2[count]=(x[i]-mx)*(x[i]-mx);
			tmp3[count]=(y[i]-my)*(y[i]-my);
			count++;
		}
	}
	tmp4=0; tmp5=0; tmp6=0;
	for(int i=0;i<count;i++)
	{
		tmp4+=tmp1[i];
		tmp5+=tmp2[i];
		tmp6+=tmp3[i];
	}
	
	
	*r=tmp4/sqrt(tmp5*tmp6);
	
	free(tmp1);
	free(tmp2);
	free(tmp3);
	
	return Asc_utils::NO_ERR;
}	

// method: scalar_sub
// 
// arguments:
//  
//

// return: a value of Asc_utils::ERR if there is an error
long Asc_utils::scalar_sub(double **input1,int rows1,int cols1,double scalar,double ***M)
{
        long error_code=0;
	
	
	    /* *M = (double **)malloc(rows1 * sizeof(double *));
		for(int i = 0; i < rows1; i++)
		{
			(*M)[i] = (double *)malloc(cols1 * sizeof(double));
		} */
		
		for (int i=0;i<rows1;i++)
			for(int j=0;j<cols1;j++)
				(*M)[i][j]=scalar-input1[i][j];
				
	 return Asc_utils::NO_ERR;	
}

 long Asc_utils::search_vec(long * vec, long len, long search_value,long ** result_index)
 {
	 //*result_index=(long*)malloc(len*sizeof(long));
	 long result_len=0;
	 for (long i=0;i<len;i++)
	 {
		if (vec[i]==search_value)
		{
			(*result_index)[result_len]=i;
			result_len++;
		}
	 }
	 
	return result_len;	 
 
 }
 
 
  long Asc_utils::sel_index_vec(long * vec, long len_vec, long *index,long index_len,long ** result)
 {
	// *result=(long*)malloc(index_len*sizeof(long));
	 
	 for (int i=0;i<index_len;i++)
	 {		
			(*result)[i]=vec[index[i]];	
	 }
	 
	return 0;	 
 
 }
 
double Asc_utils::mean_vec(double *vec,long len)
 {
	
	long count=0;
	double result=0;
	
	for(long i=0;i<len;i++)
	{
		if (isnan(vec[i])==0)
		{
			result+=vec[i];
			count++;
		}
	}
	//printf("count %ld %f ",count,result);
	result=result/(double)count;
	return result;
	
 }
 
 double Asc_utils::normal_rand()
 {
    
	double u = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double v = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double r = u * u + v * v;
    if (r == 0 || r > 1) return Asc_utils::normal_rand();
    double c = sqrt(-2 * log(r) / r);
    return (u * c);
    //return ((double)rand()/RAND_MAX);
 }
 
 // median
// simple function to compute  median and  ignore nan values
double Asc_utils::median_vec(double *arr, int n) 
{
double median;
long count=0;
for(int i=0;i<n;i++)
{
	if (isnan(arr[i])==0)
		count++;
}
n=count;
std::sort(&arr[0], &arr[n]);
 median=(n % 2 ? arr[n / 2] : (arr[n / 2 - 1] + arr[n / 2]) / 2);
 
 return median;
}

//result=input/sum(abs(input))
long Asc_utils::divide_sumabs(double *input,long len,double **result)
{
	double sumabs=0;
    //*result=(double*)malloc(len*sizeof(double));
	for(int i=0;i<len;i++)
		sumabs+=(sqrt(input[i]*input[i]));
	for(int i=0;i<len;i++)
		(*result)[i]=input[i]/sumabs; 
 return 0;
}


 
 
 
 