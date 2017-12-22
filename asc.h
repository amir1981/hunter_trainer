// file: asc.h
//

// make sure definitions are only made once
//
#ifndef ISIP_ASC
#define ISIP_ASC

// include standard libraries
//
#include <stdio.h>         // basic io stuff
#include <stdlib.h>        // basic C++ stuff
#include <sys/stat.h>      // filename processing
#include <string.h>        // C string functions
#include <time.h>


#include "asc_model.h"
#include "spectrogram/spectrogram.h"
#include "asc_utils/asc_utils.h"

// Asc: a signal processing class that performs automated sound comparison
// using a simple dynamic time-warping approach. This technology was
// originally developing in MATLAB, but then converted to C++ for
// efficiency and integration into a web application.
//

//define window types; currencly just rectangular and Hammings
#define RECT  0
#define HAMMING 1


class Asc {
  
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
  static const long ERR = 99999;
  
  //---------------------------------------------------------------------------
  //
  // protected data
  //
  //---------------------------------------------------------------------------
protected:

  // define a buffer for audio data:
  //  This is equivalent to 30 seconds x 8,000 samples per second
  //
  static const long MAX_AUDIO_SIZE = 240000;
  
  //input signal
  long num_samples;
  double* sig;
  
  //model class
  Asc_model a_model;
 
  //configuration data
  
	short nFFT; //FFT size

	short winLen; //window length

	short nOverlap; // number of overlaped samples

	short winType; // Window type
	
	short minTrain; //minimum training  files


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
  ~Asc();

  // method: default constructor
  //
  Asc();

  //---------------------------------------------------------------------------
  //
  // other public methods
  //
  //---------------------------------------------------------------------------
public:

  // computational methods:
  //  enroll is used to initialize a model; update is used
  //  to update an existing model given a new audio file;
  //  evaluate is used to compare a test file to an existing model.
  //
  int train(int nf_a, char** fnames, char* mname);
  int evaluate(int nf_a, char** fnames, char* mname);
  int info(char *mname_a);
  
  int load_config(char *file_name);

  //---------------------------------------------------------------------------
  //
  // private methods
  //
  //---------------------------------------------------------------------------
private:

  // model manipulation methods:
  //  There are two groups of functions. The first three are used
  //  for training. The second set are used for evaluation.
  //
  double init_model(char* afname);
  double update_model(char* afname);
  double  score1(double *inp_sig,long len);
  double  score2(double *inp_sig,long len);
  double score3(double *inp_sig,long len);
  double compute_score(double *inp_sig,long len);
  
  int save_model(char* mname);

  int load_model(char* name);
  float evaluate_model(char* afname);

  // load audio file:
  //  reads a sampled data file in "raw" format
  //
  int load_audio(char* af_e);

  //
  // end of class
};

// end of include file
//
#endif
