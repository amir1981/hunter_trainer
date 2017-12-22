// file: spectrogram.h
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
#include <math.h>          //math functions

// kiss FFT library 
#include "kiss_fft130/kiss_fft.h"

#define PI 3.14159265
//define window types; currencly just rectangular and Hammings
#define RECT  0
#define HAMMING 1

//define output_type options
#define NORMAL  0
#define ABS     1
#define POW     2
#define LOGPOW  3

// Spectrogram: A class that performs spectrogram operation on the input signal
//
class Spectrogram {
  
  //--------------------------------------------------------------------------
  //
  // public constants
  //
  //--------------------------------------------------------------------------
public:
  
  //where calculated spectrogram is saved
  //each column contains nfreq_bins
  //each row contains ntime_bins	

  //if  out_type is normal spectrogram1 contains real part and spectrogram2 contains imiginary part
  //in other cases spectrogram2 is NULL and spectrogram1 conatins the spectrogram  
 double **spectrogram1; //real part 
                      				  
 double **spectrogram2; //imiginary part
 
 double **spectrogram_abs; //in NORMAL case abs
 
 double *mean_select_re; //selected cols mean
 double *mean_select_im;
 
 


  //number of rows and columns of spectrogram
  long nFreq_bins; // number of frequncy bins = nFFT/2+1 (if nFFT is even)
 
  long nTime_bins; //number of time bins; calculated after moving the window over the signal
 
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
  static const long ERR_SIGPARAMs=-3; //input signal or params are not set
  static const long ERR = 9999;
  
  //---------------------------------------------------------------------------
  //
  // protected data
  //
  //---------------------------------------------------------------------------
protected:

double *input_signal; //where input signal is saved

double *seg; // segment of the input signal




//fft real and imiginary parts
double *fft_re;
double *fft_im;

long num_samples; // number of samples in the input signal (length of input vector)



short nFFT; //FFT size

short winLen; //window length

short nOverlap; // number of overlaped samples

short winType; // Window type

double *win; // window coefficients    

short outputType; //output type (NORMAL,ABS,POW,LOGPOW)


 //  used by  kiss fft lib
 //size_t buflen;
 kiss_fft_cpx  * in;
 kiss_fft_cpx  * out;
 kiss_fft_cfg  cfg;
  

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
  ~Spectrogram();

  // method: default constructor
  //
  Spectrogram();
  //Spectrogram(double *sig,long num_samples);

  //---------------------------------------------------------------------------
  //
  // other public methods
  //
  //---------------------------------------------------------------------------
public:
 
  
  //set input signal
  long set_signal(double* sig, long num_samples);
  //set parameters for spectrogram
  long set_params(short nFFT,short winLen,short nOverlap,short winType,short outputType);
 
  //calculate spectrogram
  long spectrogram(void);

  //mean of selected Cols : rsult saved under mean_select_re and *_im
  long mean_selected_cols(long *cols,long len);

  // return the spectrogram
  //if  out_type is normal spectrogram1 contains real part and spectrogram2 contains imiginary part
  //in other cases spectrogram2 is NULL and spectrogram1 conatins the spectrogram
  //int return_spectrogram(double **spectrogram1,double **spectrogram2); 

  //---------------------------------------------------------------------------
  //
  // private methods
  //
  //---------------------------------------------------------------------------
private:
	//calculate FFT
	long fft_init(void);
	long fft(double *seg); 
	//set  "win"  variable based on window type
	long set_window(void);
	//calculate and return nTime_bins
	long cal_nTimebins(void); 

  

  //
  // end of class
};

// end of include file
//
//#endif
