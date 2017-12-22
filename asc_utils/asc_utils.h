// file: asc_utls.h
//

// make sure definitions are only made once
//
//#ifndef ISIP_ASC
//#define ISIP_ASC

// include standard libraries
//
#include <stdio.h>         // basic io stuff
#include <string.h>
#include <stdlib.h>        // basic C++ stuff
#include <sys/stat.h>      // filename processing
#include <iostream>        // string functions

#include <math.h>
#include <limits.h>
#include <float.h>
#include <algorithm>


using namespace std;

// asc_utils: Contains utility functions used in asc project
// all functions are defined as static cases
class Asc_utils {
  
  //--------------------------------------------------------------------------
  //
  // public constants
  //
  //--------------------------------------------------------------------------
public:
  
  // define the class name
  //
  static const char* CLASS_NAME;

  

  //----------------------------------------
  //
  // error codes
  //
  //----------------------------------------  

  static const long NO_ERR = 0;
  static const long ERR_FILE = -1;
  static const long ERR_MODEL = -1;
  static const long ERR_MEM=-2;
  static const long ERR = 99999;
  
  //---------------------------------------------------------------------------
  //
  // protected data
  //
  //---------------------------------------------------------------------------
protected:

  // define variables to hold and parse filenames
  //
  std::string mdir;
  std::string sid;
  std::string pname;

  

  //---------------------------------------------------------------------------
  //
  // required public methods
  //
  //---------------------------------------------------------------------------
public:
  
  // method: name
  //
  static const char* name() {
    return CLASS_NAME;
  }

  // method: destructor
  //
  ~Asc_utils();

  // method: default constructor
  //
  Asc_utils();

  //---------------------------------------------------------------------------
  //
  // other public methods
  //
  //---------------------------------------------------------------------------
public:
   static void test();
  
   //Use dynamic programming to find min cost path through matrix M
   static long DP_path(long **p,long **q,double **M,int m,int n);
   
   //array_min
   static long array_min(double *input,long len);
   
   //array sqrt
   static long array_sqrt(double *input,long len,double **result);
   
   //Matrix multiply
   static long matrix_multiply(double **input1,int rows1,int cols1,double **input2,int rows2,int cols2,double ***result);
   
   //Matrix sum
   static long matrix_sum(double **input,int rows,int cols,int dim,double **result);
   
   //Matrix point pow
   static long matrix_epow2(double **input,int rows,int cols,double ***result);
   
   //set NaN members
   static long matrix_setNaN(double **input,int rows,int cols,double set_to,double *** result);
   
   //matrix copy
   static long matrix_copy(double **input,int rows,int cols,double ***result);
   
   //array copy
   static long array_copy(double *input,int len,double **result);
   
   //matrix_sim
   static long matrix_sim(double **input1,int rows1,int cols1,double **input2,int rows2,int cols2,double ***M);
   
   //subtract a scalar from matrix
   static long scalar_sub(double **input1,int rows1,int cols1,double scalar,double ***M);
   
   //correlation
   static long linear_corr(double *x,double* y,int len,double *r);
   
   //search in vecotr
   static long search_vec(long * vec, long len, long search_value,long ** result_index);

   //select indiced in a vector
   static long sel_index_vec(long * vec, long len_vec, long *index,long index_len,long ** result);
   
   //
   static double mean_vec(double *vec,long len);
   
   static double normal_rand();
   
   static double median_vec(double *arr, int n);
   
   static  long divide_sumabs(double *input,long len,double **result);
   
 static long array_pow2(double *input,long len,double ** result);
   
  //---------------------------------------------------------------------------
  //
  // private methods
  //
  //---------------------------------------------------------------------------
private:



  //
  // end of class
};

// end of include file
//
//#endif
