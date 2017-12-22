// file: asc_2.cc
//
// this file contains public methods
//

// isip include files
//
#include "asc.h"
#define _OPEN_SYS
#include <dirent.h>
#undef _OPEN_SYS

//------------------------------------------------------------------------
//
// public methods
//
//-----------------------------------------------------------------------

// method: train
// 
// arguments:
//  char** fnames: array of filenames
//  char* mname: model filename (output)
//
// return: an int flag indicating status
//
// This method creates a model from several audio files. Note that scores
// for each file are printed to stdout.
//
int Asc::train(int nf_a, char** fnames_a, char* mname_a) {

  // declare local variables
  //
  float score;
  int status;
  int model_existed=0;
 
 
  // use the first file to initialize the model
  //
  if (nf_a < 1) {
    return Asc::ERR_FILE;
  }

  //create  models  dir if not existed


  FILE *istream;
  if ( (istream = fopen ( mname_a, "r" ) ) == NULL )
  {
	model_existed=0;
  }
  else
  {
	model_existed=1;
	fclose(istream);
  }

//if the model is not existed
if (model_existed==0)
{	
    
	//initilize the model with the first file
	
     score = init_model(fnames_a[0]);
	  fprintf(stdout, "%f\n",score);
	  if (score < 0) {
	    
		return Asc::ERR_MODEL;
	  }
	 else{
	
			for (int i=1; i<nf_a; i++) //read files
			{
				score = update_model(fnames_a[i]);
				fprintf(stdout, "%f \n",score);
				if (score < 0) {
				  return Asc::ERR_MODEL;
				}
				

			}
		
			save_model(mname_a);
		
       }

}
else //model existed
{
	
	//load the model
	if (Asc::load_model(mname_a) != Asc::NO_ERR) {
    fprintf(stdout, "*> Asc::evaluate: error loading model (%s)\n",
	    mname_a);
    return Asc::ERR_FILE;
    }
	
	for (int i=0; i<nf_a; i++) //read files
	{
		score = update_model(fnames_a[i]);
		fprintf(stdout, "%f\n",score);
		if (score < 0) {
		  return Asc::ERR_MODEL;
		}
	
		
	
	}
	

			save_model(mname_a);

			
}
  
  
///////////  
   
  // exit gracefully
  //
  return Asc::NO_ERR;
}

// method: evaluate
// 
// arguments:
//  char** fnames: array of filenames
//  char* mname: model filename (output)
//
// return: an int flag indicating status
//
// This method evaluates a file against an existing model.
//
int Asc::evaluate(int nf_a, char** fnames_a, char* mname_a) {

  // declare local variables
  //
  double score;
  int status;

  // load a model
  //  
if (Asc::load_model(mname_a) != Asc::NO_ERR) {
    fprintf(stdout, "*> Asc::evaluate: error loading model (%s)\n",
	    mname_a);
    return Asc::ERR_FILE;
  }

  // loop over the audio files
  //
  for (int i=0; i<nf_a; i++) {

    // evaluate each file
    //
    score = evaluate_model(fnames_a[i]);
    fprintf(stdout, "%f\n",score);
    if (score < 0) {
      return Asc::ERR_MODEL;
    }
  }

  // exit gracefully
  //
  return Asc::NO_ERR;
}

//Show info for a model
int Asc::info(char *mname_a)
{
    // load a model
	  //
	  std::string model_dir;
	  std::string model_file;
	  
	  model_dir.append("models/");
	  model_dir.append(mname_a);
	  model_file.append(model_dir); 
	  model_file.append("/");
	  model_file.append(mname_a);
	if (Asc::load_model((char*)model_file.c_str()) != Asc::NO_ERR) {
    fprintf(stdout, "*> Asc::evaluate: error loading model (%s)\n",
	    mname_a);
    return Asc::ERR_FILE;
  }
  fprintf(stdout,"model name :%s\n",mname_a);
  fprintf(stdout,"number of records :%ld\n\n",node::len);
  
  a_model.retrieve();
  
  return Asc::NO_ERR;
}

//load from  config  file
//curently just set the variables but should be  implemented upon the requiest.
int Asc::load_config(char *file_name)
{
	//temp  later  should loaded from  file
	
	nFFT=256; //FFT size

	winLen=120; //window length

	nOverlap=40; // number of overlaped samples

	winType=HAMMING; // Window type
	
	minTrain=3; //minimum training  files
	
	 // exit gracefully
	//
	return Asc::NO_ERR;
}