// file: asc.cc
//
// modified:
//  20130610: redefined the command line interface
//

// isip include files
//
#include "asc.h"

// asc: automated sound comparison utility
//
// This is a driver program that demonstrates the functionality
// of the sound comparison module.
//
int main(int argc, const char** argv) {

  // create a sound comparison object
  //
  Asc asc;
  srand(time(NULL));
  
  // declare local variables
  //
  long status = Asc::NO_ERR;
  double score;
  

  
  srand(time(NULL));
  //---------------------------------------------------------------------------
  // configure and parse the command line:
  //  nothing fancy here, a fixed argument structure is assumed
  //---------------------------------------------------------------------------

  // check the arguments
  //
  if (argc < 3) {
    system("cat help_message.text");
    exit(-1);
  }

  // grab the arguments
  //
  char* mode = (char*)argv[1];
  char* model_name = (char*)argv[2];
  int nf = argc - 3;

  char* fnames[nf];
  int j = 0;
  for (int i=3; i<argc; i++) {
    fnames[j++] = (char*)argv[i];
  }

  //--------------------------------------------------------------------------
  //
  // main processing loop
  //
  //---------------------------------------------------------------------------

  //load config file
  //config.txt contains  configuration 
  // for  now  it contains  the  config in itself.
  asc.load_config((char *)"");
  
  // check for training
  //  
  if (strcmp(argv[1], "train") == 0) {

    // train a model
    //
    status = asc.train(nf, fnames, model_name);
  }

  // check for evaluation
  //
  else if (strcmp(argv[1], "evaluate") == 0) {

    // evaluate a model
    //
    status = asc.evaluate(nf, fnames, model_name);
  }
   else if (strcmp(argv[1], "info") == 0) {

    // evaluate a model
    //
    status = asc.info(model_name);
  }
  // else: unknown mode
  //
  else {
    status = Asc::ERR;
  }
  
  // display status
  //
  if (status != Asc::NO_ERR) {
    fprintf(stdout, "%f\n", (float)status);
    exit(Asc::ERR);
  }
  
  // exit gracefully
  //
  exit(Asc::NO_ERR);
}
