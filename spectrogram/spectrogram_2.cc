// file: spectrogram_2.cc
//
// this file contains public methods
//

// isip include files
//
#include "spectrogram.h"


//------------------------------------------------------------------------
//
// public methods
//
//-----------------------------------------------------------------------

// method: set_siganl
// 
// arguments:
//  double *sig: input signal
//  long num_samples:  number of samples of the input signal
//
// return: a value of Spectrogram::ERR_FILE if there is an error
//
// This method sets the input signal and its length
//
long Spectrogram::set_signal(double* sig, long num_samples) {

  // create signal buffer
  
  this->num_samples = num_samples;
  this->input_signal = (double *)malloc(num_samples*sizeof(double));
 
  for (int i=0;i<num_samples;i++)
	this->input_signal[i]=sig[i];
	
  // exit gracefully
  //
  return Spectrogram::NO_ERR;
}


long Spectrogram::set_params(short nFFT,short winLen,short nOverlap,short winType,short outputType)
{
	long error_code=0;
	
	if (this->num_samples == 0)
	{
		return Spectrogram::ERR_SIGPARAMs;
	}
	this->nFFT=nFFT;
	this->winLen=winLen;
	this->nOverlap=nOverlap;
	this->winType=winType;
	this->outputType=outputType;
    
	error_code=set_window();
	 
	error_code=cal_nTimebins();
	
	this->nFreq_bins=(short)(nFFT/2)+1;
    
	fft_re= (double *) malloc(nFFT*sizeof(double));
	fft_im= (double *) malloc(nFFT*sizeof(double));
	
	seg= (double *) malloc(winLen*sizeof(double));
	 
	//check for error in allocating memory
	if (fft_re==NULL || fft_im==NULL || seg==NULL)
		return Spectrogram::ERR_MEM;
		
	//allocating 2 dimensional spectrograms
	if (outputType==NORMAL) //allocate both spectrograms1 &2
	{
		spectrogram1 = (double **)malloc(nFreq_bins* sizeof(double *));
		for (int i = 0; i < nFreq_bins; i++) 
			spectrogram1[i] = (double *)malloc(nTime_bins * sizeof(double));
	    
		spectrogram2 = (double **)malloc(nFreq_bins* sizeof(double *));
		for (int i = 0; i < nFreq_bins; i++) 
			spectrogram2[i] = (double *)malloc(nTime_bins * sizeof(double));	
			
		spectrogram_abs = (double **)malloc(nFreq_bins* sizeof(double *));
		for (int i = 0; i < nFreq_bins; i++) 
			spectrogram_abs[i] = (double *)malloc(nTime_bins * sizeof(double));
			
		if (spectrogram1 ==NULL || spectrogram2==NULL || spectrogram_abs==NULL)
			return Spectrogram::ERR_MEM;	
	}
	else
	{
		spectrogram1 = (double **)malloc(nFreq_bins* sizeof(double *));
		for (int i = 0; i < nFreq_bins; i++) 
			spectrogram1[i] = (double *)malloc(nTime_bins * sizeof(double));
		if (spectrogram1 ==NULL)
			return Spectrogram::ERR_MEM;	
	}

	 mean_select_re = (double *)malloc(nFreq_bins* sizeof(double ));
		
	 mean_select_im = (double *)malloc(nFreq_bins* sizeof(double ));
		
			
    //allocate memory for fft
	error_code=fft_init();
	
	 
	if (error_code != Spectrogram::NO_ERR)
			return error_code;	

   // exit gracefully
  //
 
  return Spectrogram::NO_ERR;

}

// method: spectrogram
// 
// arguments:

// return: a value of Spectrogram::ERR if there is an error
//
//calculate the spectrogram
long Spectrogram::spectrogram(void)
{
	long index=0;
	int tbin=0;
	int k,j;
	long error_code;
	
	
	
	while (index+winLen-1<=num_samples)
	{
		for (k=0,j=index;j<=index+winLen-1;j++,k++)
			seg[k]=input_signal[j]*win[k]; // windowed signal
		error_code=fft(seg);
		
		
		for (int i=0;i<nFreq_bins;i++)
		{
			switch(outputType)
			{
				case NORMAL:
					spectrogram1[i][tbin]=fft_re[i];
					spectrogram2[i][tbin]=fft_im[i];
					spectrogram_abs[i][tbin]=sqrt(spectrogram1[i][tbin]*spectrogram1[i][tbin]+spectrogram2[i][tbin]*spectrogram2[i][tbin]);
				break;
				case ABS:
					spectrogram1[i][tbin]=sqrt(fft_re[i]*fft_re[i]+fft_im[i]*fft_im[i]);
				break;
				case POW:
					spectrogram1[i][tbin]=fft_re[i]*fft_re[i]+fft_im[i]*fft_im[i];
				break;
				case LOGPOW:
					spectrogram1[i][tbin]=log10(fft_re[i]*fft_re[i]+fft_im[i]*fft_im[i]);
				break; 
			}
		}
		index=index+winLen-nOverlap;
		tbin++;
		
	}

}
long Spectrogram::mean_selected_cols(long *cols,long len)
{
	double * sum_re,*sum_im;
	sum_re=(double*)malloc(nFreq_bins*sizeof(double));
    sum_im=(double*)malloc(nFreq_bins*sizeof(double));
	
	
	
	for (int i=0;i<nFreq_bins;i++)
	{
		sum_re[i]=0;
		sum_im[i]=0;
	}
	for (int k=0;k<len;k++)
	{
		for (int i=0;i<nFreq_bins;i++)
		{
		  if(isnan(spectrogram1[i][cols[k]])==0 && isnan(spectrogram1[i][cols[k]])==0)
		  {
			sum_re[i]=sum_re[i]+spectrogram1[i][cols[k]];	
			sum_im[i]=sum_im[i]+spectrogram2[i][cols[k]];
			
		  }
		 }
	 }
	 
	 for (int i=0;i<nFreq_bins;i++)
	 {
		mean_select_re[i]=sum_re[i]/len;
	    mean_select_im[i]=sum_im[i]/len;
		/* if(isnan(mean_select_re[i])==1)
		printf("NAN\n"); */
	 
	 }	
    
	 
	 
	free(sum_re);
    free(sum_im);
	return 0;
	 
}
