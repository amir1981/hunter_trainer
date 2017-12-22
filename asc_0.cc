// file: asc_0.cc
//

// isip include files
//
#include "asc.h"

//------------------------------------------------------------------------
//
// basic required methods
//
//-----------------------------------------------------------------------

// method: default constructor
//
Asc::Asc() {

  // create an audio buffer
  //
  num_samples = 0;
  sig = (double*)malloc(Asc::MAX_AUDIO_SIZE*sizeof(double));

  // exit gracefully
  //
};

// method: destructor
//
Asc::~Asc() {

  // clean up memory used to store the signal
  //
  free(sig);

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
const char* Asc::CLASS_NAME("Asc");
