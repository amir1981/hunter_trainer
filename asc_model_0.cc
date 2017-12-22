// file: asc_model_0.cc
//

// isip include files
//
#include "asc_model.h"

//------------------------------------------------------------------------
//
// basic required methods
//
//-----------------------------------------------------------------------

// method: default constructor
//

node::node(double **spectrogram_re,double **spectrogram_im,long nRows,long nCols,char * filename,node * next)
{
    
	this->spectrogram_re = (double **)malloc(nRows* sizeof(double *));
		for (int i = 0; i < nRows; i++) 
			this->spectrogram_re[i] = (double *)malloc(nCols * sizeof(double));

	this->spectrogram_im = (double **)malloc(nRows* sizeof(double *));
		for (int i = 0; i < nRows; i++) 
			this->spectrogram_im[i] = (double *)malloc(nCols * sizeof(double));
		
	for (long i=0;i<nRows;i++)
		for	(long j=0;j<nCols;j++)
			this->spectrogram_re[i][j]=spectrogram_re[i][j];
	
	for (long i=0;i<nRows;i++)
		for	(long j=0;j<nCols;j++)
			this->spectrogram_im[i][j]=spectrogram_im[i][j];

	
			
    strcpy(this->filename,filename);
	len++;
	example_no=len;
	
	this->next=next;

}

// method: destructor
//
node::~node()
{
	 for (int i = 0; i < nRows; i++) 
		free(spectrogram_re[i]);
     free(spectrogram_re);
	 
	  for (int i = 0; i < nRows; i++) 
		free(spectrogram_im[i]);
     free(spectrogram_im);
}

// method: default constructor
//

Asc_model::Asc_model() {
	head = NULL;
	
  // exit gracefully
  //
}

// method: destructor
//
Asc_model::~Asc_model() {

for (int i = 0; i < nRows; i++) 
		free(first_abs[i]);
    free(first_abs);
	
for (int i = 0; i < nRows; i++) 
		free(mean_re[i]);
    free(mean_re);

for (int i = 0; i < nRows; i++) 
		free(mean_im[i]);
    free(mean_im);
	
for (int i = 0; i < nRows; i++) 
		free(mean_abs[i]);
    free(mean_abs);	
	
for (int i = 0; i < nRows; i++) 
		free(abs_mean[i]);
    free(abs_mean);	
	
for (int i = 0; i < nRows; i++) 
		free(mean_timebin_re[i]);
    free(mean_timebin_re);	

for (int i = 0; i < nRows; i++) 
		free(mean_timebin_im[i]);
    free(mean_timebin_im);	
	
for (int i = 0; i < nRows; i++) 
		free(mean_freqbin_re[i]);
    free(mean_freqbin_re);
	
for (int i = 0; i < nRows; i++) 
		free(mean_freqbin_im[i]);
    free(mean_freqbin_im);	

	
	
node* current = head;
while( current != 0 ) {
    node* next = current->next;
    delete current;
    current = next;
}
head = 0;
  
  //
}

//-----------------------------------------------------------------------------
//
// we define non-integral constants in the default constructor
//      
//-----------------------------------------------------------------------------

// constants: class name
//
const char* Asc_model::CLASS_NAME("Asc_model");

//initilize static members
long node::len=0;
long node::nRows=0;
long node::nCols=0;
long Asc_model::nRows=0;
long Asc_model::nCols=0;
double **Asc_model::first_abs=NULL;
double **Asc_model::mean_re=NULL;
double **Asc_model::mean_im=NULL;
double **Asc_model::mean_abs=NULL;
double **Asc_model::abs_mean=NULL;
double **Asc_model::mean_timebin_re=NULL;
double **Asc_model::mean_timebin_im=NULL;
double **Asc_model::mean_freqbin_re=NULL;
double **Asc_model::mean_freqbin_im=NULL;
double Asc_model::A1_score1=-100;
double Asc_model::B1_score1=-100;
double Asc_model::A1_score2=-100;
double Asc_model::B1_score2=-100;
double Asc_model::A1_score3=-100;
double Asc_model::B1_score3=-100;

