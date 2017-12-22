// file: spectrogram.cc
//

// isip include files
//
#include "spectrogram.h"

//------------------------------------------------------------------------
//
// basic required methods
//
//-----------------------------------------------------------------------

// method: default constructor
//
Spectrogram::Spectrogram() {

  // 
  //
  this->num_samples = 0;
  this->input_signal = NULL;
  this->nFFT=0;
  this->winLen=0; 
  this->nOverlap=0;
  this->winType=0;
  this->nFreq_bins=0;
  this->nTime_bins=0;
  
  fft_re=NULL;
  fft_im=NULL;
  spectrogram1=NULL;
  spectrogram2=NULL;
  spectrogram_abs=NULL;
  mean_select_re=NULL;
  mean_select_im=NULL;
  
  seg=NULL;
  
  in=NULL;
  out=NULL;
  cfg=NULL;
  
  // exit gracefully
  //
};

/* Spectrogram::Spectrogram(double *sig,long num_samples) {

  // create signal buffer
  //
  this->num_samples = num_samples;
  this->input_signal = (double *)malloc(num_samples*sizeof(double));
  for (int i=0;i<num_samples;i++)
	this->input_signal[i]=sig[i];

  fft_re=NULL;
  fft_im=NULL;
  spectrogram_re=NULL;
  spectrogram_im=NULL;
  spectrogram_abs=NULL;
  spectrogram_pow=NULL;
  spectrogram_logpow=NULL;
  
  // exit gracefully
  //
}; */



// method: destructor
//
Spectrogram::~Spectrogram() {

  // clean up memory used to store the signal
  //
  free(input_signal); //free  input_signal
  free(win); //free win
  free(seg);
  
  free(fft_re);
  free(fft_im);
	
 free(mean_select_re);
  
 
  free(mean_select_im);	
  //free 2-d arrays	

 
  
  for (int i = 0; i < nFreq_bins; i++) 
	free(spectrogram1[i]);
  free(spectrogram1);
  
  if(outputType == NORMAL)
  {
	for (int i = 0; i < nFreq_bins; i++) 
		free(spectrogram2[i]);
	free(spectrogram2);
	
	 for (int i = 0; i < nFreq_bins; i++) 
	free(spectrogram_abs[i]);
  free(spectrogram_abs);
  }
  
    // free kiss fft stuff here  
    free(in);
    free(out);
    free(cfg);

  // exit gracefully
  //
};

//-----------------------------------------------------------------------------
//
// we define non-integral constants in the default constructor
//      
//-----------------------------------------------------------------------------

// constants: class name
//
const char* Spectrogram::CLASS_NAME("Spectrogram");


