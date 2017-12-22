// file: spectrogram_test.cc
//

// isip include files
//
#include "spectrogram.h"

// spectrogram: automated sound comparison utility
//
// This is a driver program that demonstrates the functionality
// of the sound comparison module.
//
long readfile(double **sig,char *filename);
int savefile(double *sig,long num_samples, char *filename);

int main(int argc, const char** argv) {

  // create a sound comparison object
  //
  Spectrogram spg;
  double *sig;
  long num_samples;
  long error_code;
  
  // declare local variables
  //
 

  //---------------------------------------------------------------------------
  // configure and parse the command line:
  //  nothing fancy here, a fixed argument structure is assumed
  //---------------------------------------------------------------------------

  // check the arguments
  //
  if (argc < 1) {
    system("cat help_message.text");
    exit(-1);
  }

  // grab the arguments
  //
  char* input_file = (char*)argv[1];

  //char* input_type = (char*)argv[3]; //short , int, double
  
 
  // echo the command line
  //
  fprintf(stdout, "This program is to test the spectrogram.\n");
  fprintf(stdout, "The input signal should be saved in \"double\" format.\n");
  fprintf(stdout, "input file: %s\n", input_file);
 
  
  //--------------------------------------------------------------------------
  //
  // main processing loop
  //
  //---------------------------------------------------------------------------
   
   num_samples=readfile(&sig,input_file);
   
   if (spg.set_signal(sig,num_samples) !=0)
   {
	fprintf(stdout," Error in set_signal.\n");
	 exit(-1);
   } 
   if ((error_code=spg.set_params(64,4,1,HAMMING,NORMAL)) !=0)
   {
	fprintf(stdout," Error in set_params %ld \n",error_code);
	 exit(-1);
   }
   spg.spectrogram();
   for (int i=0;i<spg.nFreq_bins;i++)
   {
	for (int j=0;j<spg.nTime_bins;j++)
	{
		fprintf(stdout, "(%f,%f)\t",spg.spectrogram1[i][j],spg.spectrogram2[i][j]);
	}
	fprintf(stdout,"\n");
   }
   
   //savefile(spg.sigpow,num_samples,output_file);

  // exit gracefully
  //
  
  free(sig);
  exit(0);
}
long readfile(double **sig,char *filename)
{
  // open the file read only
  //
  long num_samples;
  FILE* fp = fopen(filename, "r");

  if (fp == NULL) {
    fprintf(stdout, "*> Spectrogram_test:: readfile: error opening file (%s)\n",
	    filename);
    return 0;
  }

  // grab the file size
  //
  fseek (fp, 0, SEEK_END);
  long num_samples_in_file = ftell(fp) / sizeof(double);
  rewind (fp);

  // load the file into memory
  //
  *sig = (double*)malloc(num_samples_in_file*sizeof(double));
  num_samples = fread(*sig, sizeof(double), num_samples_in_file, fp);
 if (num_samples != num_samples_in_file) {
 fprintf(stdout, "*> Spectrogram_test:: readfile: Number of samples is not the same as expected.");
 return 0;
 }
  // close the file
  //
  fclose(fp);

  return num_samples;
}
int savefile(double *sig,long num_samples, char *filename)
{
	FILE *fo=fopen(filename,"w");
	for(int i=0;i<num_samples;i++)
	{
		fprintf(fo,"%f ",sig[i]);
	}	
	fclose(fo);	
	return 0;
}


