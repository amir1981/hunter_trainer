// file: spectrogram_1.cc
//
// This file contains implementations of the private methods.
//


// isip include files
//
#include "spectrogram.h"



//------------------------------------------------------------------------
//
// private methods
//
//-----------------------------------------------------------------------

// method: fft_init
// 
// arguments:

// return: a value of Spectrogram::ERR if there is an error
//
// Initilize and alloacte fft memory
long Spectrogram::fft_init(void)
{
 
	
    in = (kiss_fft_cpx *) malloc(nFFT *sizeof(kiss_fft_cpx));
    out= (kiss_fft_cpx *) malloc(nFFT *sizeof(kiss_fft_cpx));
    cfg = kiss_fft_alloc(nFFT,0,0,0);
	
	if (in == NULL || out== NULL || cfg == NULL)
		return Spectrogram::ERR_MEM;
		
	// exit gracefully
    //
    return Spectrogram::NO_ERR;	
}


// method: fft
// 
// arguments:

// return: a value of Spectrogram::ERR if there is an error
//
// Calculate the fft of input_signal

//This function used "kiss fft"  lib but it is possible to use other  libraries.
long Spectrogram::fft(double * seg)
{
    if (in == NULL)
		return Spectrogram::ERR;
	
    //zero padding	
	for (int k=0;k<nFFT;++k) {
        in[k].r = 0;
        in[k].i = 0;
    }
	for (int k=0;k<winLen;++k) {
        in[k].r = seg[k];
        in[k].i = 0;
    }
	//call kiss fft
	kiss_fft(cfg,in,out);
	
	
	for (int k=0;k<nFFT;++k) {
		fft_re[k]=out[k].r;
		fft_im[k]=out[k].i;
    }
	
	
    // exit gracefully
    //
    return Spectrogram::NO_ERR;
}

// method: set_window
// 
// arguments:

// return: a value of Spectrogram::ERR if there is an error
//
// Calculate the window coefficients and allocate the memory 

long Spectrogram::set_window(void)
{
	//allocate the memory for the window
	win=(double *)malloc(winLen*sizeof(double));
	if (win==NULL) return Spectrogram::ERR_MEM;
	
	if (winType == RECT)
	{
		for (int i=0;i<winLen;i++)
			win[i]=1;  //rect window
			
		return Spectrogram::NO_ERR;
	}
	if (winType == HAMMING)
	{
		for (int i=0;i<winLen;i++)
			win[i]=0.54-0.46*cos(2*PI*i/(winLen-1)); //hamming window
	
		return Spectrogram::NO_ERR;
	}
	//else
	//return Spectrogram::ERR; 
}

// method: cal_nTimebins
// 
// arguments:

// return: a value of Spectrogram::ERR if there is an error
//
//calculate and return nTime_bins
long Spectrogram::cal_nTimebins(void)
{
	nTime_bins=0;
	int i=0;
	while(i+winLen-1<=num_samples)
	{	
		i=i+winLen-nOverlap;
		nTime_bins++;
	}
	return Spectrogram::NO_ERR;
}

