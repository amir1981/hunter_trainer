// file: asc_1.cc
//
// This file contains implementations of the private methods.
//

// isip include files
//
#include "asc.h"

//------------------------------------------------------------------------
//
// private methods
//
//-----------------------------------------------------------------------

// method: init_model
// 
// arguments:
//  char* afname_a: input audio filename
//
// return: a score indicating the quality of the model
//
// This method initializes a model using a sample audio file.
//
double Asc::init_model(char* afname_a) {

  // declare local variables
  //
  long status = Asc::NO_ERR;
  double score = 100.0; // The score for the first example is always 100 at initiliztion step


  //define spectrogram
  Spectrogram spg;
  
   
   
  // load the audio
  //
  if ((status = Asc::load_audio(afname_a)) != Asc::NO_ERR) {
    fprintf(stdout, "*> Asc:: init_model: error loading file (%s)\n",
	    afname_a);
    //return(-1000);
	exit(-1);
  }    
  //set signal   
  if (spg.set_signal(sig,num_samples) !=0)
   {
	fprintf(stdout," Error in set_signal.\n");
	return Asc::ERR;
   } 
   //set params
   if ((status=spg.set_params(nFFT,winLen,nOverlap,HAMMING,NORMAL)) !=0)
   {
	fprintf(stdout," Error in set_params %ld \n",status);
	 return Asc::ERR;
   }
   
   //calculate the spectrogram
   spg.spectrogram();
   
   //initilize the model
   a_model.init(spg.nFreq_bins,spg.nTime_bins);
  
   //add example
   a_model.add_example(spg.spectrogram1,spg.spectrogram2,spg.nFreq_bins,spg.nTime_bins,afname_a);
   
   a_model.abs_first();
  // exit gracefully
  //
  return score;
}

// method: update_model
// 
// arguments:
//  char* afname: input audio file
//
// return: a score indicating the quality of the model
//
// This method implements an algorithm to update a model.
//
double Asc::update_model(char* afname_a) {

  // declare local variables
  //
  long status;
  double score = 0;
  long mn;
  double **B_new_re=NULL;
  double **B_new_im=NULL;
  
  
     B_new_re = (double **)malloc((Asc_model::nRows)* sizeof(double *));
	   for (int i = 0; i < (Asc_model::nRows); i++) 
			B_new_re[i] = (double *)malloc( (Asc_model::nCols) * sizeof(double));
	  B_new_im = (double **)malloc((Asc_model::nRows)* sizeof(double *));
	    for (int i = 0; i < (Asc_model::nRows); i++) 
		    B_new_im[i] = (double *)malloc( (Asc_model::nCols) * sizeof(double));
			
			
//

   Asc_model::A1_score1=-100;
 Asc_model::B1_score1=-100;
  Asc_model::A1_score2=-100;
 Asc_model::B1_score2=-100;
  Asc_model::A1_score3=-100;
 Asc_model::B1_score3=-100;

 
  //define spectrogram
  Spectrogram spg;
  
  
  // load the audio
  //
  if ((status = Asc::load_audio(afname_a)) != Asc::NO_ERR) {
    fprintf(stdout, "*> Asc:: update_model: error loading file (%s)\n",
	    afname_a);
    //return Asc::ERR_FILE;
	exit(-1);
  }    

 //compute the score before presenting the new training file
 score=compute_score(sig,num_samples);


 
	  //set signal   
	  if (spg.set_signal(sig,num_samples) !=0)
	   {
		fprintf(stdout," Error in set_signal.\n");
		return Asc::ERR;
	   } 
	   //set params
	   if ((status=spg.set_params(nFFT,winLen,nOverlap,HAMMING,NORMAL)) !=0)
	   {
		fprintf(stdout," Error in set_params %ld \n",status);
		 return Asc::ERR;
	   }
	   
	   //calculate the spectrogram
	   spg.spectrogram();
	  
	   
	   {
			   double **SM=NULL;
			   double **M=NULL;
			   long *p;
			   long *q;
			   
				SM = (double **)malloc(Asc_model::nCols * sizeof(double *));
				for(int i = 0; i < Asc_model::nCols; i++)
					SM[i] = (double *)malloc(spg.nTime_bins * sizeof(double));
					
				 M = (double **)malloc(Asc_model::nCols * sizeof(double *));
				for(int i = 0; i < Asc_model::nCols; i++)
					M[i] = (double *)malloc(spg.nTime_bins * sizeof(double));	
					
				p = (long *)malloc(100*Asc_model::nCols*spg.nTime_bins* sizeof(long));
				q = (long *)malloc(100*Asc_model::nCols*spg.nTime_bins* sizeof(long));
				
				
			   if ((status=Asc_utils::matrix_sim(Asc_model::first_abs,Asc_model::nRows,Asc_model::nCols,spg.spectrogram_abs,spg.nFreq_bins,spg.nTime_bins,&SM)) !=0)
			   {
				 fprintf(stdout," Error in matrix_sim %ld \n",status);
				 return Asc::ERR;
			   
			   }		
			   
					
				Asc_utils::scalar_sub(SM,Asc_model::nCols,spg.nTime_bins,1,&M);
				
				
				
				if ((status=Asc_utils::DP_path(&p,&q,M,Asc_model::nCols,spg.nTime_bins)) <0)
			   {
				 fprintf(stdout," Error in DP_path %ld \n",status);
				 return Asc::ERR;
			   
			   }
			   mn=status; 
			  
			   
			   		
				{
					long *index=NULL;
					long *qr=NULL;
					index=(long*)malloc(mn*sizeof(long));
					qr=(long*)malloc(mn*sizeof(long)); // allocate more memory (instead of index_len, mn)
					
					for (long j=0;j<Asc_model::nCols;j++)
					{		
						long index_len=Asc_utils::search_vec(p,mn,j,&index);
						
						Asc_utils::sel_index_vec(q,mn,index,index_len,&qr);
						
						spg.mean_selected_cols(qr,index_len);
						
						for(int i=0;i<Asc_model::nRows;i++)
						{
							B_new_re[i][j]=spg.mean_select_re[i];
							B_new_im[i][j]=spg.mean_select_im[i];
						}
						
					}
					free(index);
					free(qr);
				}
				  //free  SM
				 for (int i = 0; i < Asc_model::nCols; i++) 
					 free(SM[i]);
				  free(SM);	
				   for (int i = 0; i < Asc_model::nCols;i++) 
					 free(M[i]);
				  free(M);
				  
				  free(p);
				  free(q); 
		}
	  // exit gracefully
	   //add example
	
	  //  if (score >0 || node::len>minTrain)
       // { 
			 a_model.add_example(B_new_re,B_new_im,Asc_model::nRows,Asc_model::nCols,afname_a);
        //}	

   
   for (int i = 0; i < Asc_model::nRows; i++) 
	 free(B_new_re[i]);
  free(B_new_re);
   for (int i = 0; i < Asc_model::nRows; i++) 
	 free(B_new_im[i]);
  free(B_new_im);
  
  return score;
}

// method: save_model
// 
// arguments:
//  char* mname_a: full filename for model file
//
// return: a flag indicating status
//
// This function writes a model to a model directory using
// a fixed filename format.
//
int Asc::save_model(char* mname_a) {


  a_model.save(mname_a);
  

  // exit gracefully
  //
  return Asc::NO_ERR;
}

// method: load_model
// 
// arguments:
//  char* mname_a: full filename of input model
//
// return: a boolean value indicating status
//
// This method reads a model from a model file. The filename
// follows the same format used in create_model().
//
int Asc::load_model(char* mname_a) {

  a_model.load(mname_a);
  a_model.mean_all();
  a_model.abs_first();
  
  // exit gracefully
 
  return Asc::NO_ERR;
}

// method: evaluate_model
// 
// arguments:
//  char* afname: input audio file
//
// return: a floating point score in the range [0, 100]. if
//         less than zero, an error has occurred.
//
// This method implements the algorithm that compares a file
// to a model and produces a score.
//
float Asc::evaluate_model(char* afname_a) {

  // declare local variables
  //
  long status = Asc::NO_ERR;
  double score = 0;

  // load the audio data
  //
  if ((status = Asc::load_audio(afname_a)) != Asc::NO_ERR) {
    fprintf(stdout, "*> Asc:: evaluate_model: error loading audio (%s)\n",
	    afname_a);
    return Asc::ERR_FILE;
  }    

  // evaluation here
  //
  
  
  //un comment if you want to re calculate the coefficients for evaluation (not  needed)
  
/*   Asc_model::A1_score1=-100;
 Asc_model::B1_score1=-100;
  Asc_model::A1_score2=-100;
 Asc_model::B1_score2=-100;
  Asc_model::A1_score3=-100;
 Asc_model::B1_score3=-100; */
  score = compute_score(sig,num_samples);
  //a_model.print();
  // exit gracefully
  //
  return score;
}

// method: load_audio
// 
// arguments:
//  char* afname: audio filename to be used for enrollment
//
// return: a value of Asc::ERR_FILE if there is an error
//
// This method loads and audio file into memory and stores it as
// protected data in the class.
//
int Asc::load_audio(char* afname_a) {

  // open the file read only
  //
  short * sig_tmp;
  sig_tmp = (short int*)malloc(Asc::MAX_AUDIO_SIZE*sizeof(short int));
  double *sig_tmp2;
  sig_tmp2 = (double*)malloc(Asc::MAX_AUDIO_SIZE*sizeof(double));
  
  
  FILE* fp = fopen(afname_a, "r");

  if (fp == NULL) {
    fprintf(stdout, "*> Asc:: load_audio: error opening file (%s)\n",
	    afname_a);
    return Asc::ERR_FILE;
  }

  // grab the file size
  //
  fseek (fp, 0, SEEK_END);
  long num_samples_in_file = ftell(fp) / sizeof(short int);
  rewind (fp);

  if (num_samples_in_file > Asc::MAX_AUDIO_SIZE) {
    fclose(fp);
    return Asc::ERR_FILE;
  }

  // load the file into memory
  //
  num_samples = fread(sig_tmp, sizeof(short int), num_samples_in_file, fp);
 
  for(int i=0;i<num_samples;i++) //convert to double  between -1 and 1
  {
	sig_tmp2[i]=(double)sig_tmp[i]/32768.0;
  }
  //divide by sum of absoulte values 
  Asc_utils::divide_sumabs(sig_tmp2,num_samples,&sig);
  
  // close the file
  //
  fclose(fp);
  free(sig_tmp);
  free(sig_tmp2);
 
  // check for consistency
  //
  if (num_samples != num_samples_in_file) {
    return Asc::ERR_FILE;
  }

  // exit gracefully
  //
  return Asc::NO_ERR;
}

double Asc::compute_score(double *inp_sig,long len)
{
 double score_1,score_2,score_3;
 double score;
 
 a_model.mean_all();
 
 score_1=score1(inp_sig,len);
 score_2=score2(inp_sig,len);
 score_3=score3(inp_sig,len);
 double tmp1[2];
 double tmp2[2];
 
 
 tmp1[0]=score_1;
 tmp1[1]=score_2;
 tmp2[0]=score_3;
 tmp2[1]=Asc_utils::mean_vec(tmp1,2);
 score=floor(100*Asc_utils::mean_vec(tmp2,2));
 return (score);
 
}
//score  using correlation  over freq  
//used  mean to  find the average   correlation
//This function compute the score and calibrate it. The calibration is done
//using the ideal reference and noisy signal. 
//calibration fourmula:
//final_score=(score-B)/(A-B);
//A=average score for ideal cases
//B=average score for the noise

double Asc::score1(double *inp_sig,long len)
{
	
	// should move and allocate outside :speedup
	double **sp_re;
	double **sp_im;
	double **sp_abs;
	
	
	//can go outside: speedup
	double *tmp1;
	double *tmp2;
	
	//speedup
	double *match_freq0;
	
	
	double *ideal_raws;
	double *noise_raws; 
	
	double *tmp3;
	double *tmp4;
	
	double raw_score=0;
	double score=0;
	
	//speedup
	long no_noise_examples=3;//number of noise examples
	
	tmp1=(double*)malloc(Asc_model::nRows*sizeof(double));
	tmp2=(double*)malloc(Asc_model::nRows*sizeof(double));
	match_freq0=(double*)malloc(Asc_model::nCols*sizeof(double));
	
	ideal_raws=(double*)malloc(node::len*sizeof(double));
	noise_raws=(double*)malloc(no_noise_examples*sizeof(double)); // no_noise_examples example of noise
	
	tmp3=(double*)malloc(len*sizeof(double));
	tmp4=(double*)malloc(len*sizeof(double));
	
	
	sp_re = (double **)malloc(Asc_model::nRows* sizeof(double *));
	for (int i = 0; i < Asc_model::nRows; i++) 
		sp_re[i] = (double *)malloc(Asc_model::nCols * sizeof(double));
	sp_im = (double **)malloc(Asc_model::nRows* sizeof(double *));
	for (int i = 0; i < Asc_model::nRows; i++) 
		sp_im[i] = (double *)malloc(Asc_model::nCols * sizeof(double));	
	sp_abs = (double **)malloc(Asc_model::nRows* sizeof(double *));
	for (int i = 0; i < Asc_model::nRows; i++) 
		sp_abs[i] = (double *)malloc(Asc_model::nCols * sizeof(double));
	
	//Calculate raw score for  ideal examples
	
 if (Asc_model::A1_score1==-100 || Asc_model::B1_score1==-100)
 {	
	for (int itr0=0;itr0<node::len;itr0++)
	{	
		a_model.return_node(itr0,&sp_re,&sp_im,&sp_abs); //reterive example itr0
			
		for (int i=0;i<Asc_model::nCols;i++)
		{
		 double r0;
			for (int j=0;j<Asc_model::nRows;j++)
			{
				tmp1[j]=sp_abs[j][i];
				tmp2[j]=Asc_model::mean_abs[j][i];
			}
		Asc_utils::linear_corr(tmp1,tmp2,Asc_model::nRows,&r0);	
	    match_freq0[i]=sqrt(r0*r0);		
		
		}
	ideal_raws[itr0]=Asc_utils::mean_vec(match_freq0,Asc_model::nCols);
	//printf("score1 ideals:%f %f\n",ideal_raws[itr0],Asc_utils::mean_vec(match_freq0,Asc_model::nCols));
	}


	
	//Calculate raw score for  noise examples
	
	// no_noise_examples example of random signal
	{
	
	      //allocation can be done outside of loop :  speedup
	      //double **SM=NULL;
		  long mn;
		 // double **M;
		  
		  
		  //move out: speedup
		  double **B_new_re;
		  double **B_new_im;
		  
		  long status;
		  
		  //move ou: speedup
		  double *match_freq1=NULL;
		  
		  match_freq1=(double*)malloc(len*sizeof(double));
		  
	       B_new_re = (double **)malloc((Asc_model::nRows)* sizeof(double *));
		   for (int i = 0; i < (Asc_model::nRows); i++) 
				B_new_re[i] = (double *)malloc( (Asc_model::nCols) * sizeof(double));
		   B_new_im = (double **)malloc((Asc_model::nRows)* sizeof(double *));
		   for (int i = 0; i < (Asc_model::nRows); i++) 
				B_new_im[i] = (double *)malloc( (Asc_model::nCols) * sizeof(double)); 
				
				
			
			
			for (int itr0=0;itr0<no_noise_examples;itr0++)
			{			
				for (int j=0;j<len;j++)
					tmp4[j]=((double)rand()/RAND_MAX);
					
				Asc_utils::divide_sumabs(tmp4,len,&tmp3); // save the tmp3  divided by its abs sum	
			
				  //define spectrogram
				  Spectrogram spg;
				 
				  //set signal   
				  if (spg.set_signal(tmp3,len) !=0)
				   {
					fprintf(stdout," Error in set_signal.\n");
					return Asc::ERR;
				   } 
				   //set params
				   if ((status=spg.set_params(nFFT,winLen,nOverlap,HAMMING,NORMAL)) !=0)
				   {
					fprintf(stdout," Error in set_params %ld \n",status);
					 return Asc::ERR;
				   }
				   
				   //calculate the spectrogram
				   spg.spectrogram();
				   
				 
				  { // We have to define it here to prevent from memory over-writing 
					  double **SM;
					  double **M;
					  long *p1=NULL;
				      long *q1=NULL;
					  p1=(long*)malloc(Asc_model::nCols*spg.nTime_bins*sizeof(long));
					  q1=(long*)malloc(Asc_model::nCols*spg.nTime_bins*sizeof(long));
				      SM = (double **)malloc(Asc_model::nCols * sizeof(double *));
					  for(int i = 0; i < Asc_model::nCols; i++)
						SM[i] = (double *)malloc(spg.nTime_bins * sizeof(double));
					  M = (double **)malloc(Asc_model::nCols * sizeof(double *));
					  for(int i = 0; i < Asc_model::nCols; i++)
						M[i] = (double *)malloc(spg.nTime_bins * sizeof(double));
					
					   if ((status=Asc_utils::matrix_sim(Asc_model::abs_mean,Asc_model::nRows,Asc_model::nCols,spg.spectrogram_abs,spg.nFreq_bins,spg.nTime_bins,&SM)) !=0)
					   {
						 fprintf(stdout," Error in matrix_sim %ld \n",status);
						 return Asc::ERR;
					   
					   } 
															   
						Asc_utils::scalar_sub(SM,Asc_model::nCols,spg.nTime_bins,1,&M);
					
						
						if ((status=Asc_utils::DP_path(&p1,&q1,M,Asc_model::nCols,spg.nTime_bins)) <0)
					    {
						 fprintf(stdout," Error in DP_path %ld \n",status);
						 return Asc::ERR;
					   
					    } 
					
					    mn=status; 
						
					  { //block 2
							long *index=NULL;
							long *qr=NULL;
							
							index=(long*)malloc(mn*sizeof(long));
							qr=(long*)malloc(mn*sizeof(long));
							
								
							for (long j=0;j<Asc_model::nCols;j++)
							{
								 
								
								
								 long index_len=Asc_utils::search_vec(p1,mn,j,&index);
								
								 Asc_utils::sel_index_vec(q1,mn,index,index_len,&qr);
								
								spg.mean_selected_cols(qr,index_len);
								
								  for(int i=0;i<Asc_model::nRows;i++)
								  {
									 B_new_re[i][j]=spg.mean_select_re[i];
									 B_new_im[i][j]=spg.mean_select_im[i];
								  }
									 
								 
							}	
							
							free(index);
							free(qr);
                       } 
				        
						free(p1);
						free(q1);
						
						for (int i = 0; i < Asc_model::nCols; i++) 
							free(SM[i]);
						free(SM);
						for (int i = 0; i < Asc_model::nCols; i++) 
						    free(M[i]);
						free(M);
						
					} //end of block 1
						
					
					
					for (int i=0;i<Asc_model::nCols;i++)
					{
					 double r0;
						for (int j=0;j<Asc_model::nRows;j++)
						{
							tmp1[j]=sqrt(B_new_re[j][i]*B_new_re[j][i]+B_new_im[j][i]*B_new_im[j][i]);
							
							tmp2[j]=Asc_model::mean_abs[j][i];
						}
						
					Asc_utils::linear_corr(tmp1,tmp2,Asc_model::nRows,&r0);	
					match_freq1[i]=sqrt(r0*r0);	
                  					
					}
				noise_raws[itr0]=Asc_utils::mean_vec(match_freq1,Asc_model::nCols);  //raw score for noise
					
			//printf("socre1: noise %f %f\n",noise_raws[itr0],Asc_utils::mean_vec(match_freq1,Asc_model::nCols));
			}
		
		 
		  
	    for (int i = 0; i < Asc_model::nRows; i++) 
		 free(B_new_re[i]);
	  free(B_new_re);
	  
	   for (int i = 0; i < Asc_model::nRows; i++) 
		 free(B_new_im[i]);
	  free(B_new_im); 
	
	 
	  
	  free(match_freq1);
	
	}
	
	   //coeficients to map the  raw score 
		
		Asc_model::A1_score1=Asc_utils::mean_vec(ideal_raws,node::len);
		Asc_model::B1_score1=Asc_utils::mean_vec(noise_raws,no_noise_examples);
	
	}
	
	//Calculate raw score for the input signal
	
	{
		  //speed up
	      double **SM=NULL;  
		  long mn;
		  double **M;
		  
		  //speedup
		  double **B_new_re=NULL;
		  double **B_new_im=NULL;
		  long status;
		  
		  //speed up
		  double *match_freq1;	  
		  match_freq1=(double*)malloc(len*sizeof(double));
		  
	      B_new_re = (double **)malloc((Asc_model::nRows)* sizeof(double *));
		   for (int i = 0; i < (Asc_model::nRows); i++) 
				B_new_re[i] = (double *)malloc( (Asc_model::nCols) * sizeof(double));
		   B_new_im = (double **)malloc((Asc_model::nRows)* sizeof(double *));
		   for (int i = 0; i < (Asc_model::nRows); i++) 
				B_new_im[i] = (double *)malloc( (Asc_model::nCols) * sizeof(double));
			
           		
				
					
				  //define spectrogram
				  Spectrogram spg;
				 
				  //set signal   
				  if (spg.set_signal(inp_sig,len) !=0)
				   {
					fprintf(stdout," Error in set_signal.\n");
					return Asc::ERR;
				   } 
				   //set params
				   if ((status=spg.set_params(nFFT,winLen,nOverlap,HAMMING,NORMAL)) !=0)
				   {
					fprintf(stdout," Error in set_params %ld \n",status);
					 return Asc::ERR;
				   }
				   
				   //calculate the spectrogram
				   spg.spectrogram();
				   {
				   
							long *p1=NULL;
							long *q1=NULL;
							p1=(long*)malloc(Asc_model::nCols*spg.nTime_bins*sizeof(long));
							q1=(long*)malloc(Asc_model::nCols*spg.nTime_bins*sizeof(long));
							
						   SM = (double **)malloc(Asc_model::nCols * sizeof(double *));
						   for(int i = 0; i < Asc_model::nCols; i++)
							 SM[i] = (double *)malloc(spg.nTime_bins * sizeof(double));
						   M = (double **)malloc(Asc_model::nCols * sizeof(double *));
						   for(int i = 0; i < Asc_model::nCols; i++)
							  M[i] = (double *)malloc(spg.nTime_bins * sizeof(double));	
						   
						   if ((status=Asc_utils::matrix_sim(Asc_model::abs_mean,Asc_model::nRows,Asc_model::nCols,spg.spectrogram_abs,spg.nFreq_bins,spg.nTime_bins,&SM)) !=0)
						   {
							 fprintf(stdout," Error in matrix_sim %ld \n",status);
							 return Asc::ERR;
						   
						   }
							Asc_utils::scalar_sub(SM,Asc_model::nCols,spg.nTime_bins,1,&M);
							
					  
					  
									  
							if ((status=Asc_utils::DP_path(&p1,&q1,M,Asc_model::nCols,spg.nTime_bins)) <0)
						   {
							 fprintf(stdout," Error in DP_path %ld \n",status);
							 return Asc::ERR;
						   
						   }
						   mn=status; 
						  
							{	
								long *index=NULL;
								long *qr=NULL;
								index=(long*)malloc(mn*sizeof(long));
								qr=(long*)malloc(mn*sizeof(long));
								
								 for (long j=0;j<Asc_model::nCols;j++)
								 {
									
									 long index_len=Asc_utils::search_vec(p1,mn,j,&index);
									
									 Asc_utils::sel_index_vec(q1,mn,index,index_len,&qr);
									
									spg.mean_selected_cols(qr,index_len);
									
									  for(int i=0;i<Asc_model::nRows;i++)
									  {
										 B_new_re[i][j]=spg.mean_select_re[i];
										 B_new_im[i][j]=spg.mean_select_im[i];
									  }
									
									
								}
								 
									 
									 free(index);
									 free(qr);
							}
							for (int i = 0; i < Asc_model::nCols; i++) 
								free(SM[i]);
							free(SM);
							for (int i = 0; i < Asc_model::nCols; i++) 
								 free(M[i]);
							free(M);		
							free(p1);
							
							free(q1);
						
					}
					
					for (int i=0;i<Asc_model::nCols;i++)
					{
					 double r0;
						for (int j=0;j<Asc_model::nRows;j++)
						{
							tmp1[j]=sqrt(B_new_re[j][i]*B_new_re[j][i]+B_new_im[j][i]*B_new_im[j][i]);
							tmp2[j]=Asc_model::mean_abs[j][i];
						}
					
						 
					Asc_utils::linear_corr(tmp1,tmp2,Asc_model::nRows,&r0);	
					match_freq1[i]=sqrt(r0*r0);		
					
					
					
					}
				raw_score=Asc_utils::mean_vec(match_freq1,Asc_model::nCols);
				//printf("score1 raw: %f %f\n",raw_score,Asc_utils::mean_vec(match_freq1,Asc_model::nCols));
					
			  

       
			
		
		
	    
	
	   for (int i = 0; i < Asc_model::nRows; i++) 
		 free(B_new_re[i]);
	  free(B_new_re);
	  
	   for (int i = 0; i < Asc_model::nRows; i++) 
		 free(B_new_im[i]);
	  free(B_new_im);
	 
	   
	  
	  
	  free(match_freq1);
	
	}
	
	 // Calculate the final score: Use average raw scores for  noise and ideal to map the raw score for signal into interval [0,1]			
			
		
		
		//mapping
		
		if (Asc_model::B1_score1 != Asc_model::A1_score1)
			score=(raw_score-Asc_model::B1_score1)/(Asc_model::A1_score1-Asc_model::B1_score1);
		else
			score=raw_score;
		
		
		
		if (score<0)
			score=0;
		if (score>1)
			score=1;
			
			 

	 for (int i = 0; i < Asc_model::nRows; i++) 
		free(sp_re[i]);
    free(sp_re);	
	 for (int i = 0; i < Asc_model::nRows; i++) 
		free(sp_im[i]);
    free(sp_im);
	 for (int i = 0; i < Asc_model::nRows; i++) 
		free(sp_abs[i]);
    free(sp_abs);
	
	free(tmp1);
	free(tmp2);
	free(match_freq0);
	free(ideal_raws);
	free(noise_raws);
	
	free(tmp3);
	free(tmp4);
  
return(score);
}
//score  using correlation  over freq  
//used  median to  find the average   correlation
//This function compute the score and calibrate it. The calibration is done
//using the ideal reference and noisy signal. 
//calibration fourmula:
//final_score=(score-B)/(A-B);
//A=average score for ideal cases
//B=average score for the noise

double Asc::score2(double *inp_sig,long len)
{
	
	// should move and allocate outside :speedup
	double **sp_re;
	double **sp_im;
	double **sp_abs;
	
	
	//can go outside: speedup
	double *tmp1;
	double *tmp2;
	
	//speedup
	double *match_freq0;
	
	
	double *ideal_raws;
	double *noise_raws; 
	
	double *tmp3;
	double *tmp4;
	
	double raw_score=0;
	double score=0;
	
	//speedup
	long no_noise_examples=3;//number of noise examples
	
	tmp1=(double*)malloc(Asc_model::nRows*sizeof(double));
	tmp2=(double*)malloc(Asc_model::nRows*sizeof(double));
	match_freq0=(double*)malloc(Asc_model::nCols*sizeof(double));
	
	ideal_raws=(double*)malloc(node::len*sizeof(double));
	noise_raws=(double*)malloc(no_noise_examples*sizeof(double)); // no_noise_examples example of noise
	
	tmp3=(double*)malloc(len*sizeof(double));
	tmp4=(double*)malloc(len*sizeof(double));
	
	
	sp_re = (double **)malloc(Asc_model::nRows* sizeof(double *));
	for (int i = 0; i < Asc_model::nRows; i++) 
		sp_re[i] = (double *)malloc(Asc_model::nCols * sizeof(double));
	sp_im = (double **)malloc(Asc_model::nRows* sizeof(double *));
	for (int i = 0; i < Asc_model::nRows; i++) 
		sp_im[i] = (double *)malloc(Asc_model::nCols * sizeof(double));	
	sp_abs = (double **)malloc(Asc_model::nRows* sizeof(double *));
	for (int i = 0; i < Asc_model::nRows; i++) 
		sp_abs[i] = (double *)malloc(Asc_model::nCols * sizeof(double));
	
	//Calculate raw score for  ideal examples
	
if (Asc_model::A1_score2==-100  ||  Asc_model::B1_score2==-100)
{

	for (int itr0=0;itr0<node::len;itr0++)
	{	
		a_model.return_node(itr0,&sp_re,&sp_im,&sp_abs); //reterive example itr0
			
		for (int i=0;i<Asc_model::nCols;i++)
		{
		 double r0;
			for (int j=0;j<Asc_model::nRows;j++)
			{
				tmp1[j]=sp_abs[j][i];
				tmp2[j]=Asc_model::mean_abs[j][i];
			}
		Asc_utils::linear_corr(tmp1,tmp2,Asc_model::nRows,&r0);	
	    match_freq0[i]=sqrt(r0*r0);		
		
		}
	ideal_raws[itr0]=Asc_utils::median_vec(match_freq0,Asc_model::nCols);
	
	}


	
	//Calculate raw score for  noise examples
	
	// no_noise_examples example of random signal
	{
	
	      //allocation can be done outside of loop :  speedup
	      //double **SM=NULL;
		  long mn;
		 // double **M;
		  
		  
		  //move out: speedup
		  double **B_new_re;
		  double **B_new_im;
		  
		  long status;
		  
		  //move ou: speedup
		  double *match_freq1=NULL;
		  
		  match_freq1=(double*)malloc(len*sizeof(double));
		  
	       B_new_re = (double **)malloc((Asc_model::nRows)* sizeof(double *));
		   for (int i = 0; i < (Asc_model::nRows); i++) 
				B_new_re[i] = (double *)malloc( (Asc_model::nCols) * sizeof(double));
		   B_new_im = (double **)malloc((Asc_model::nRows)* sizeof(double *));
		   for (int i = 0; i < (Asc_model::nRows); i++) 
				B_new_im[i] = (double *)malloc( (Asc_model::nCols) * sizeof(double)); 
				
				
			
			
			for (int itr0=0;itr0<no_noise_examples;itr0++)
			{			
				for (int j=0;j<len;j++)
					tmp4[j]=((double)rand()/RAND_MAX);
					
				Asc_utils::divide_sumabs(tmp4,len,&tmp3); // save the tmp3  divided by its abs sum	
			
				  //define spectrogram
				  Spectrogram spg;
				 
				  //set signal   
				  if (spg.set_signal(tmp3,len) !=0)
				   {
					fprintf(stdout," Error in set_signal.\n");
					return Asc::ERR;
				   } 
				   //set params
				   if ((status=spg.set_params(nFFT,winLen,nOverlap,HAMMING,NORMAL)) !=0)
				   {
					fprintf(stdout," Error in set_params %ld \n",status);
					 return Asc::ERR;
				   }
				   
				   //calculate the spectrogram
				   spg.spectrogram();
				   
				 
				  { // We have to define it here to prevent from memory over-writing 
					  double **SM;
					  double **M;
					  long *p1=NULL;
				      long *q1=NULL;
					  p1=(long*)malloc(Asc_model::nCols*spg.nTime_bins*sizeof(long));
					  q1=(long*)malloc(Asc_model::nCols*spg.nTime_bins*sizeof(long));
				      SM = (double **)malloc(Asc_model::nCols * sizeof(double *));
					  for(int i = 0; i < Asc_model::nCols; i++)
						SM[i] = (double *)malloc(spg.nTime_bins * sizeof(double));
					  M = (double **)malloc(Asc_model::nCols * sizeof(double *));
					  for(int i = 0; i < Asc_model::nCols; i++)
						M[i] = (double *)malloc(spg.nTime_bins * sizeof(double));
					
					   if ((status=Asc_utils::matrix_sim(Asc_model::abs_mean,Asc_model::nRows,Asc_model::nCols,spg.spectrogram_abs,spg.nFreq_bins,spg.nTime_bins,&SM)) !=0)
					   {
						 fprintf(stdout," Error in matrix_sim %ld \n",status);
						 return Asc::ERR;
					   
					   } 
															   
						Asc_utils::scalar_sub(SM,Asc_model::nCols,spg.nTime_bins,1,&M);
					
						
						if ((status=Asc_utils::DP_path(&p1,&q1,M,Asc_model::nCols,spg.nTime_bins)) <0)
					    {
						 fprintf(stdout," Error in DP_path %ld \n",status);
						 return Asc::ERR;
					   
					    } 
					
					    mn=status; 
						
					  { //block 2
							long *index=NULL;
							long *qr=NULL;
							
							index=(long*)malloc(mn*sizeof(long));
							qr=(long*)malloc(mn*sizeof(long));
							
								
							for (long j=0;j<Asc_model::nCols;j++)
							{
								 
								
								
								 long index_len=Asc_utils::search_vec(p1,mn,j,&index);
								
								 Asc_utils::sel_index_vec(q1,mn,index,index_len,&qr);
								
								spg.mean_selected_cols(qr,index_len);
								
								  for(int i=0;i<Asc_model::nRows;i++)
								  {
									 B_new_re[i][j]=spg.mean_select_re[i];
									 B_new_im[i][j]=spg.mean_select_im[i];
								  }
									 
								 
							}	
							
							free(index);
							free(qr);
                       } 
				        
						free(p1);
						free(q1);
						
						for (int i = 0; i < Asc_model::nCols; i++) 
							free(SM[i]);
						free(SM);
						for (int i = 0; i < Asc_model::nCols; i++) 
						    free(M[i]);
						free(M);
						
					} //end of block 1
						
					
					
					for (int i=0;i<Asc_model::nCols;i++)
					{
					 double r0;
						for (int j=0;j<Asc_model::nRows;j++)
						{
							tmp1[j]=sqrt(B_new_re[j][i]*B_new_re[j][i]+B_new_im[j][i]*B_new_im[j][i]);
							
							tmp2[j]=Asc_model::mean_abs[j][i];
						}
						
					Asc_utils::linear_corr(tmp1,tmp2,Asc_model::nRows,&r0);	
					match_freq1[i]=sqrt(r0*r0);	
                  					
					}
				noise_raws[itr0]=Asc_utils::median_vec(match_freq1,Asc_model::nCols);  //raw score for noise
					
			
			}
		
		 
		  
	    for (int i = 0; i < Asc_model::nRows; i++) 
		 free(B_new_re[i]);
	  free(B_new_re);
	  
	   for (int i = 0; i < Asc_model::nRows; i++) 
		 free(B_new_im[i]);
	  free(B_new_im); 
	
	 
	  
	  free(match_freq1);
	
	}
	
	//coeficients to map the  raw score 
	
		Asc_model::A1_score2=Asc_utils::mean_vec(ideal_raws,node::len);
		Asc_model::B1_score2=Asc_utils::mean_vec(noise_raws,no_noise_examples);
	
}	
	
	//Calculate raw score for the input signal
	
	{
		  //speed up
	      double **SM=NULL;  
		  long mn;
		  double **M;
		  
		  //speedup
		  double **B_new_re=NULL;
		  double **B_new_im=NULL;
		  long status;
		  
		  //speed up
		  double *match_freq1;	  
		  match_freq1=(double*)malloc(len*sizeof(double));
		  
	      B_new_re = (double **)malloc((Asc_model::nRows)* sizeof(double *));
		   for (int i = 0; i < (Asc_model::nRows); i++) 
				B_new_re[i] = (double *)malloc( (Asc_model::nCols) * sizeof(double));
		   B_new_im = (double **)malloc((Asc_model::nRows)* sizeof(double *));
		   for (int i = 0; i < (Asc_model::nRows); i++) 
				B_new_im[i] = (double *)malloc( (Asc_model::nCols) * sizeof(double));
			
           		
				
					
				  //define spectrogram
				  Spectrogram spg;
				 
				  //set signal   
				  if (spg.set_signal(inp_sig,len) !=0)
				   {
					fprintf(stdout," Error in set_signal.\n");
					return Asc::ERR;
				   } 
				   //set params
				   if ((status=spg.set_params(nFFT,winLen,nOverlap,HAMMING,NORMAL)) !=0)
				   {
					fprintf(stdout," Error in set_params %ld \n",status);
					 return Asc::ERR;
				   }
				   
				   //calculate the spectrogram
				   spg.spectrogram();
				   {
				   
							long *p1=NULL;
							long *q1=NULL;
							p1=(long*)malloc(Asc_model::nCols*spg.nTime_bins*sizeof(long));
							q1=(long*)malloc(Asc_model::nCols*spg.nTime_bins*sizeof(long));
							
						   SM = (double **)malloc(Asc_model::nCols * sizeof(double *));
						   for(int i = 0; i < Asc_model::nCols; i++)
							 SM[i] = (double *)malloc(spg.nTime_bins * sizeof(double));
						   M = (double **)malloc(Asc_model::nCols * sizeof(double *));
						   for(int i = 0; i < Asc_model::nCols; i++)
							  M[i] = (double *)malloc(spg.nTime_bins * sizeof(double));	
						   
						   if ((status=Asc_utils::matrix_sim(Asc_model::abs_mean,Asc_model::nRows,Asc_model::nCols,spg.spectrogram_abs,spg.nFreq_bins,spg.nTime_bins,&SM)) !=0)
						   {
							 fprintf(stdout," Error in matrix_sim %ld \n",status);
							 return Asc::ERR;
						   
						   }
							Asc_utils::scalar_sub(SM,Asc_model::nCols,spg.nTime_bins,1,&M);
							
					  
					  
									  
							if ((status=Asc_utils::DP_path(&p1,&q1,M,Asc_model::nCols,spg.nTime_bins)) <0)
						   {
							 fprintf(stdout," Error in DP_path %ld \n",status);
							 return Asc::ERR;
						   
						   }
						   mn=status; 
						  
							{	
								long *index=NULL;
								long *qr=NULL;
								index=(long*)malloc(mn*sizeof(long));
								qr=(long*)malloc(mn*sizeof(long));
								
								 for (long j=0;j<Asc_model::nCols;j++)
								 {
									
									 long index_len=Asc_utils::search_vec(p1,mn,j,&index);
									
									 Asc_utils::sel_index_vec(q1,mn,index,index_len,&qr);
									
									spg.mean_selected_cols(qr,index_len);
									
									  for(int i=0;i<Asc_model::nRows;i++)
									  {
										 B_new_re[i][j]=spg.mean_select_re[i];
										 B_new_im[i][j]=spg.mean_select_im[i];
									  }
									
									
								}
								 
									 
									 free(index);
									 free(qr);
							}
							for (int i = 0; i < Asc_model::nCols; i++) 
								free(SM[i]);
							free(SM);
							for (int i = 0; i < Asc_model::nCols; i++) 
								 free(M[i]);
							free(M);		
							free(p1);
							
							free(q1);
						
					}
					
					for (int i=0;i<Asc_model::nCols;i++)
					{
					 double r0;
						for (int j=0;j<Asc_model::nRows;j++)
						{
							tmp1[j]=sqrt(B_new_re[j][i]*B_new_re[j][i]+B_new_im[j][i]*B_new_im[j][i]);
							tmp2[j]=Asc_model::mean_abs[j][i];
						}
					
						 
					Asc_utils::linear_corr(tmp1,tmp2,Asc_model::nRows,&r0);	
					match_freq1[i]=sqrt(r0*r0);		
					
					
					
					}
				raw_score=Asc_utils::median_vec(match_freq1,Asc_model::nCols);
					
			  

       
			
		
		
	    
	
	   for (int i = 0; i < Asc_model::nRows; i++) 
		 free(B_new_re[i]);
	  free(B_new_re);
	  
	   for (int i = 0; i < Asc_model::nRows; i++) 
		 free(B_new_im[i]);
	  free(B_new_im);
	 
	   
	  
	  
	  free(match_freq1);
	
	}
	
	 // Calculate the final score: Use average raw scores for  noise and ideal to map the raw score for signal into interval [0,1]			
			
		
		
		//mapping
		
		if (Asc_model::B1_score2 != Asc_model::A1_score2)
			score=(raw_score-Asc_model::B1_score2 )/(Asc_model::A1_score2 -Asc_model::B1_score2 );
		else
			score=raw_score;
		
		
		
		if (score<0)
			score=0;
		if (score>1)
			score=1;
			
			 

	 for (int i = 0; i < Asc_model::nRows; i++) 
		free(sp_re[i]);
    free(sp_re);	
	 for (int i = 0; i < Asc_model::nRows; i++) 
		free(sp_im[i]);
    free(sp_im);
	 for (int i = 0; i < Asc_model::nRows; i++) 
		free(sp_abs[i]);
    free(sp_abs);
	
	free(tmp1);
	free(tmp2);
	free(match_freq0);
	free(ideal_raws);
	free(noise_raws);
	
	free(tmp3);
	free(tmp4);
  
return(score);
}


//DR. Obeid Moified Score
//This score is based on MSE but it uses the automatic calibration instead
//of polynomial fitting
//calibration fourmula:
//final_score=(score-B)/(A-B);
//A=average score for ideal cases
//B=average score for the noise

double Asc::score3(double *inp_sig,long len)
{
	
	// should move and allocate outside :speedup
	double **sp_re;
	double **sp_im;
	double **sp_abs;
	
	
	//can go outside: speedup
	double *tmp1;
	double *tmp2;
	
	//speedup
	double *match;
	
	
	double *ideal_raws;
	double *noise_raws; 
	
	double *tmp3;
	double *tmp4;
	
	double raw_score=0;
	double score=0;
	
	//speedup
	long no_noise_examples=3;//number of noise examples
	
	tmp1=(double*)malloc(Asc_model::nRows*sizeof(double));
	tmp2=(double*)malloc(Asc_model::nRows*sizeof(double));
	match=(double*)malloc(Asc_model::nCols*sizeof(double));
	
	ideal_raws=(double*)malloc(node::len*sizeof(double));
	noise_raws=(double*)malloc(no_noise_examples*sizeof(double)); // no_noise_examples example of noise
	
	tmp3=(double*)malloc(len*sizeof(double));
	tmp4=(double*)malloc(len*sizeof(double));
	
	
	sp_re = (double **)malloc(Asc_model::nRows* sizeof(double *));
	for (int i = 0; i < Asc_model::nRows; i++) 
		sp_re[i] = (double *)malloc(Asc_model::nCols * sizeof(double));
	sp_im = (double **)malloc(Asc_model::nRows* sizeof(double *));
	for (int i = 0; i < Asc_model::nRows; i++) 
		sp_im[i] = (double *)malloc(Asc_model::nCols * sizeof(double));	
	sp_abs = (double **)malloc(Asc_model::nRows* sizeof(double *));
	for (int i = 0; i < Asc_model::nRows; i++) 
		sp_abs[i] = (double *)malloc(Asc_model::nCols * sizeof(double));
	
	//Calculate raw score for  ideal examples
	

if (Asc_model::A1_score3==-100  || Asc_model::B1_score3==-100)
{
	

	for (int itr0=0;itr0<node::len;itr0++)
	{	
		a_model.return_node(itr0,&sp_re,&sp_im,&sp_abs); //reterive example itr0
			
		for (int i=0;i<Asc_model::nCols;i++)
		{
		 double r0;
		 for (int j=0;j<Asc_model::nRows;j++)
			{
				tmp1[j]=sp_abs[j][i]-Asc_model::mean_abs[j][i];
			}
			
			Asc_utils::array_pow2(tmp1,Asc_model::nRows,&tmp2);
			r0=Asc_utils::mean_vec(tmp2,Asc_model::nRows);
			
		
	    match[i]=sqrt(r0);	
		
		
		
		}
	ideal_raws[itr0]=100-10000*Asc_utils::mean_vec(match,Asc_model::nCols);
	//printf("%f %f\n",ideal_raws[itr0],Asc_utils::mean_vec(match,Asc_model::nCols));
	}


	
	//Calculate raw score for  noise examples
	
	// no_noise_examples example of random signal
	{
	
	      //allocation can be done outside of loop :  speedup
	      //double **SM=NULL;
		  long mn;
		 // double **M;
		  
		  
		  //move out: speedup
		  double **B_new_re;
		  double **B_new_im;
		  
		  long status;
		  
		 
		  
	       B_new_re = (double **)malloc((Asc_model::nRows)* sizeof(double *));
		   for (int i = 0; i < (Asc_model::nRows); i++) 
				B_new_re[i] = (double *)malloc( (Asc_model::nCols) * sizeof(double));
		   B_new_im = (double **)malloc((Asc_model::nRows)* sizeof(double *));
		   for (int i = 0; i < (Asc_model::nRows); i++) 
				B_new_im[i] = (double *)malloc( (Asc_model::nCols) * sizeof(double)); 
				
				
			
			
			for (int itr0=0;itr0<no_noise_examples;itr0++)
			{			
				for (int j=0;j<len;j++)
					tmp4[j]=((double)rand()/RAND_MAX);
					
				Asc_utils::divide_sumabs(tmp4,len,&tmp3); // save the tmp3  divided by its abs sum	
			
				  //define spectrogram
				  Spectrogram spg;
				 
				  //set signal   
				  if (spg.set_signal(tmp3,len) !=0)
				   {
					fprintf(stdout," Error in set_signal.\n");
					return Asc::ERR;
				   } 
				   //set params
				   if ((status=spg.set_params(nFFT,winLen,nOverlap,HAMMING,NORMAL)) !=0)
				   {
					fprintf(stdout," Error in set_params %ld \n",status);
					 return Asc::ERR;
				   }
				   
				   //calculate the spectrogram
				   spg.spectrogram();
				   
				 
				  { // We have to define it here to prevent from memory over-writing 
					  double **SM;
					  double **M;
					  long *p1=NULL;
				      long *q1=NULL;
					  p1=(long*)malloc(Asc_model::nCols*spg.nTime_bins*sizeof(long));
					  q1=(long*)malloc(Asc_model::nCols*spg.nTime_bins*sizeof(long));
				      SM = (double **)malloc(Asc_model::nCols * sizeof(double *));
					  for(int i = 0; i < Asc_model::nCols; i++)
						SM[i] = (double *)malloc(spg.nTime_bins * sizeof(double));
					  M = (double **)malloc(Asc_model::nCols * sizeof(double *));
					  for(int i = 0; i < Asc_model::nCols; i++)
						M[i] = (double *)malloc(spg.nTime_bins * sizeof(double));
					
					   if ((status=Asc_utils::matrix_sim(Asc_model::abs_mean,Asc_model::nRows,Asc_model::nCols,spg.spectrogram_abs,spg.nFreq_bins,spg.nTime_bins,&SM)) !=0)
					   {
						 fprintf(stdout," Error in matrix_sim %ld \n",status);
						 return Asc::ERR;
					   
					   } 
															   
						Asc_utils::scalar_sub(SM,Asc_model::nCols,spg.nTime_bins,1,&M);
					
						
						if ((status=Asc_utils::DP_path(&p1,&q1,M,Asc_model::nCols,spg.nTime_bins)) <0)
					    {
						 fprintf(stdout," Error in DP_path %ld \n",status);
						 return Asc::ERR;
					   
					    } 
					
					    mn=status; 
						
					  { //block 2
							long *index=NULL;
							long *qr=NULL;
							
							index=(long*)malloc(mn*sizeof(long));
							qr=(long*)malloc(mn*sizeof(long));
							
								
							for (long j=0;j<Asc_model::nCols;j++)
							{
								 
								
								
								 long index_len=Asc_utils::search_vec(p1,mn,j,&index);
								
								 Asc_utils::sel_index_vec(q1,mn,index,index_len,&qr);
								
								spg.mean_selected_cols(qr,index_len);
								
								  for(int i=0;i<Asc_model::nRows;i++)
								  {
									 B_new_re[i][j]=spg.mean_select_re[i];
									 B_new_im[i][j]=spg.mean_select_im[i];
								  }
									 
								 
							}	
							
							free(index);
							free(qr);
                       } 
				        
						free(p1);
						free(q1);
						
						for (int i = 0; i < Asc_model::nCols; i++) 
							free(SM[i]);
						free(SM);
						for (int i = 0; i < Asc_model::nCols; i++) 
						    free(M[i]);
						free(M);
						
					} //end of block 1
						
					
					
					for (int i=0;i<Asc_model::nCols;i++)
					{
					 double r0;
					 for (int j=0;j<Asc_model::nRows;j++)
						{
							tmp1[j]=sqrt(B_new_re[j][i]*B_new_re[j][i]+B_new_im[j][i]*B_new_im[j][i])-Asc_model::mean_abs[j][i];
						}
						
					Asc_utils::array_pow2(tmp1,Asc_model::nRows,&tmp2);
			        r0=Asc_utils::mean_vec(tmp2,Asc_model::nRows);
					
					match[i]=sqrt(r0);	
					
                  					
					}
				noise_raws[itr0]=100-10000*Asc_utils::mean_vec(match,Asc_model::nCols);  //raw score for noise
					//printf("%f %f\n",noise_raws[itr0],Asc_utils::mean_vec(match,Asc_model::nCols));
			
			}
		
		 
		  
	    for (int i = 0; i < Asc_model::nRows; i++) 
		 free(B_new_re[i]);
	  free(B_new_re);
	  
	   for (int i = 0; i < Asc_model::nRows; i++) 
		 free(B_new_im[i]);
	  free(B_new_im); 
	
	 
	  
	 
	
	}
	
	//coeficients to map the  raw score 
		Asc_model::A1_score3=Asc_utils::mean_vec(ideal_raws,node::len);
		Asc_model::B1_score3=Asc_utils::mean_vec(noise_raws,no_noise_examples);
	
}	
	
	
	//Calculate raw score for the input signal
	
	{
		  //speed up
	      double **SM=NULL;  
		  long mn;
		  double **M;
		  
		  //speedup
		  double **B_new_re=NULL;
		  double **B_new_im=NULL;
		  long status;
		  
		  
	      B_new_re = (double **)malloc((Asc_model::nRows)* sizeof(double *));
		   for (int i = 0; i < (Asc_model::nRows); i++) 
				B_new_re[i] = (double *)malloc( (Asc_model::nCols) * sizeof(double));
		   B_new_im = (double **)malloc((Asc_model::nRows)* sizeof(double *));
		   for (int i = 0; i < (Asc_model::nRows); i++) 
				B_new_im[i] = (double *)malloc( (Asc_model::nCols) * sizeof(double));
			
           		
				
					
				  //define spectrogram
				  Spectrogram spg;
				 
				  //set signal   
				  if (spg.set_signal(inp_sig,len) !=0)
				   {
					fprintf(stdout," Error in set_signal.\n");
					return Asc::ERR;
				   } 
				   //set params
				   if ((status=spg.set_params(nFFT,winLen,nOverlap,HAMMING,NORMAL)) !=0)
				   {
					fprintf(stdout," Error in set_params %ld \n",status);
					 return Asc::ERR;
				   }
				   
				   //calculate the spectrogram
				   spg.spectrogram();
				   {
				   
							long *p1=NULL;
							long *q1=NULL;
							p1=(long*)malloc(Asc_model::nCols*spg.nTime_bins*sizeof(long));
							q1=(long*)malloc(Asc_model::nCols*spg.nTime_bins*sizeof(long));
							
						   SM = (double **)malloc(Asc_model::nCols * sizeof(double *));
						   for(int i = 0; i < Asc_model::nCols; i++)
							 SM[i] = (double *)malloc(spg.nTime_bins * sizeof(double));
						   M = (double **)malloc(Asc_model::nCols * sizeof(double *));
						   for(int i = 0; i < Asc_model::nCols; i++)
							  M[i] = (double *)malloc(spg.nTime_bins * sizeof(double));	
						   
						   if ((status=Asc_utils::matrix_sim(Asc_model::abs_mean,Asc_model::nRows,Asc_model::nCols,spg.spectrogram_abs,spg.nFreq_bins,spg.nTime_bins,&SM)) !=0)
						   {
							 fprintf(stdout," Error in matrix_sim %ld \n",status);
							 return Asc::ERR;
						   
						   }
							Asc_utils::scalar_sub(SM,Asc_model::nCols,spg.nTime_bins,1,&M);
							
					  
					  
									  
							if ((status=Asc_utils::DP_path(&p1,&q1,M,Asc_model::nCols,spg.nTime_bins)) <0)
						   {
							 fprintf(stdout," Error in DP_path %ld \n",status);
							 return Asc::ERR;
						   
						   }
						   mn=status; 
						  
							{	
								long *index=NULL;
								long *qr=NULL;
								index=(long*)malloc(mn*sizeof(long));
								qr=(long*)malloc(mn*sizeof(long));
								
								 for (long j=0;j<Asc_model::nCols;j++)
								 {
									
									 long index_len=Asc_utils::search_vec(p1,mn,j,&index);
									
									 Asc_utils::sel_index_vec(q1,mn,index,index_len,&qr);
									
									spg.mean_selected_cols(qr,index_len);
									
									  for(int i=0;i<Asc_model::nRows;i++)
									  {
										 B_new_re[i][j]=spg.mean_select_re[i];
										 B_new_im[i][j]=spg.mean_select_im[i];
									  }
									
									
								}
								 
									 
									 free(index);
									 free(qr);
							}
							for (int i = 0; i < Asc_model::nCols; i++) 
								free(SM[i]);
							free(SM);
							for (int i = 0; i < Asc_model::nCols; i++) 
								 free(M[i]);
							free(M);		
							free(p1);
							
							free(q1);
						
					}
					
					for (int i=0;i<Asc_model::nCols;i++)
					{
					 double r0;
					 for (int j=0;j<Asc_model::nRows;j++)
						{
							tmp1[j]=sqrt(B_new_re[j][i]*B_new_re[j][i]+B_new_im[j][i]*B_new_im[j][i])-Asc_model::mean_abs[j][i];
						}
					
					Asc_utils::array_pow2(tmp1,Asc_model::nRows,&tmp2);
			        r0=Asc_utils::mean_vec(tmp2,Asc_model::nRows);
					match[i]=sqrt(r0);	
					
					
					}
					
				raw_score=100-10000*Asc_utils::mean_vec(match,Asc_model::nCols);
					//printf("%f %f\n",raw_score,Asc_utils::mean_vec(match,Asc_model::nCols));

	
	   for (int i = 0; i < Asc_model::nRows; i++) 
		 free(B_new_re[i]);
	  free(B_new_re);
	  
	   for (int i = 0; i < Asc_model::nRows; i++) 
		 free(B_new_im[i]);
	  free(B_new_im);
	 
	
	}
	
	 // Calculate the final score: Use average raw scores for  noise and ideal to map the raw score for signal into interval [0,1]			
			
		
		
		//mapping
		
		if (Asc_model::B1_score3 != Asc_model::A1_score3)
			score=(raw_score-Asc_model::B1_score3)/(Asc_model::A1_score3-Asc_model::B1_score3);
		else
			score=raw_score;
		
		
		
		if (score<0)
			score=0;
		if (score>1)
			score=1;
			
			 

	 for (int i = 0; i < Asc_model::nRows; i++) 
		free(sp_re[i]);
    free(sp_re);	
	 for (int i = 0; i < Asc_model::nRows; i++) 
		free(sp_im[i]);
    free(sp_im);
	 for (int i = 0; i < Asc_model::nRows; i++) 
		free(sp_abs[i]);
    free(sp_abs);
	
	free(tmp1);
	free(tmp2);
	free(match);
	free(ideal_raws);
	free(noise_raws);
	
	free(tmp3);
	free(tmp4);
  
return(score);
}