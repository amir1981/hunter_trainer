// file: asc_model.h
//

// make sure definitions are only made once
//
//#ifndef ISIP_ASC
//#define ISIP_ASC

// include standard libraries
//
#include <stdio.h>         // basic io stuff
#include <stdlib.h>        // basic C++ stuff
#include <sys/stat.h>      // filename processing
#include <string.h>        // C string functions
#include <math.h>

// Asc_model: Model class is a linked list to stroe all data for a model and related functions

//

class node{
public:
    node()
	{
	}
	
	node(double **spectrogram_re,double **spectrogram_im,long nRows,long nCols,char * filename,node * next);
	~node();
	
	static long set_size(long nRows,long nCols);
	
	
	node *next;
	
	
public:
  
  double **spectrogram_re;
  double **spectrogram_im;
  //size of spectrograms
  static long nRows;
  static long nCols;
  static long len; //total number of examples
  char filename[1000]; // raw file for this data
  long example_no; // no of this example
  
  //sum of all examples
 /*  static double **sum_re;
  static double **sum_im;
  static long **non_nan_len_re;
  statilong **non_nan_len_im; */
  
  
};


class Asc_model {
  
  //--------------------------------------------------------------------------
  //
  // public constants
  //
  //--------------------------------------------------------------------------
public:
  
  // define the class name
  //
  
 static double **first_abs;

 
 static double **mean_re;
 static double **mean_im;
 
 static double **mean_abs;
 static double **abs_mean;
 
 static double **mean_timebin_re;
 static double **mean_timebin_im;
  
 static double **mean_freqbin_re;
 static double **mean_freqbin_im;
 
 static long nRows;
 static long nCols;
 static const char* CLASS_NAME;
 //used for  calibration
 static double A1_score1;
 static double B1_score1;
 static double A1_score2;
 static double B1_score2;
 static double A1_score3;
 static double B1_score3;

  //----------------------------------------
  //
  // error codes
  //
  //----------------------------------------  

  static const long NO_ERR = 0;
  static const long ERR_FILE = -1;
  static const long ERR_MODEL = -1;
  static const long ERR = 99999;
  
  //---------------------------------------------------------------------------
  //
  // protected data
  //
  //---------------------------------------------------------------------------
protected:
 node *head;

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
  ~Asc_model();
  
  //Initilize static varibales
  void init(long nRows,long nCols);
  
  // method: default constructor
  //
  Asc_model();

  //---------------------------------------------------------------------------
  //
  // other public methods
  //
  //---------------------------------------------------------------------------
public:

  //add new example
  long add_example(double **spectrogram_re,double **spectrogram_im,long nRows,long nCols,char * filename);
  //mean of all examples
  long mean_all();  
  // return abs of first example
  long abs_first();
  
  long return_node(long node_no,double ***spectrogram_re,double ***spectrogram_im,double ***spectrogram_abs);
  
  long save(char *filename);
  long load(char *mname_a);
  
  long print();// print  model matrices into  screen (used for model inspection and debugging)
  long retrieve();
  
  
  //---------------------------------------------------------------------------
  //
  // private methods
  //
  //---------------------------------------------------------------------------
private:

  
};

// end of include file
//
//#endif
