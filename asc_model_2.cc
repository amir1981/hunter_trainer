// file: asc_model_2.cc
//
// this file contains public methods
//

// isip include files
//
#include "asc_model.h"

//------------------------------------------------------------------------
//
// public methods
//
//-----------------------------------------------------------------------

long node::set_size(long nR,long nC)
{	
	nRows=nR;
	nCols=nC;
	len=0;
	return 0;
}

//return a node
long Asc_model::return_node(long node_no,double ***spectrogram_re,double ***spectrogram_im,double ***spectrogram_abs)
{
	node *p;
	p=head;
	long count=0;
	
	
	
	/* 
	*spectrogram_re = (double **)malloc(Asc_model::nRows* sizeof(double *));
	for (int i = 0; i < Asc_model::nRows; i++) 
		(*spectrogram_re)[i] = (double *)malloc(Asc_model::nCols * sizeof(double));
	*spectrogram_im = (double **)malloc(Asc_model::nRows* sizeof(double *));
	for (int i = 0; i < Asc_model::nRows; i++) 
		(*spectrogram_im)[i] = (double *)malloc(Asc_model::nCols * sizeof(double));	
	*spectrogram_abs = (double **)malloc(Asc_model::nRows* sizeof(double *));
	for (int i = 0; i < Asc_model::nRows; i++) 
		(*spectrogram_abs)[i] = (double *)malloc(Asc_model::nCols * sizeof(double));
		*/	
	
if (node_no>node::len-1)
	{
		fprintf(stdout," node_no should be in the range of 0 < %ld \n",node::len-1);
		return Asc_model::ERR;
	}
	
	while(count<node_no)
	{
		p=p->next;
		count++;
	}
	for (long i=0;i<Asc_model::nRows;i++)
		for (long j=0;j<Asc_model::nCols;j++)
		{
			(*spectrogram_re)[i][j]=p->spectrogram_re[i][j];
			(*spectrogram_im)[i][j]=p->spectrogram_im[i][j];
			(*spectrogram_abs)[i][j]=sqrt(p->spectrogram_re[i][j]*p->spectrogram_re[i][j]+p->spectrogram_im[i][j]*p->spectrogram_im[i][j]);
		}
}

//add examples
long Asc_model::add_example(double **spectrogram_re,double **spectrogram_im,long nRows,long nCols,char * filename)
{
	node *p;
	if (head ==NULL)
		head= new node(spectrogram_re,spectrogram_im,nRows,nCols,filename,NULL);
	else
	{
		p=head;
		while (p->next != NULL)
			p=p->next;
		p->next=new node(spectrogram_re,spectrogram_im,nRows,nCols,filename,NULL);	
	}
	
	return Asc_model::NO_ERR;
}

//initilize static  stuff

void Asc_model::init(long nRows,long nCols)
{

    first_abs = (double **)malloc(nRows* sizeof(double *));
	for (int i = 0; i < nRows; i++) 
		first_abs[i] = (double *)malloc(nCols * sizeof(double));
		
	mean_re = (double **)malloc(nRows* sizeof(double *));
	for (int i = 0; i < nRows; i++) 
		mean_re[i] = (double *)malloc(nCols * sizeof(double));
		
	mean_im = (double **)malloc(nRows* sizeof(double *));
	for (int i = 0; i < nRows; i++) 
		mean_im[i] = (double *)malloc(nCols * sizeof(double));	
		
	mean_abs = (double **)malloc(nRows* sizeof(double *));
	for (int i = 0; i < nRows; i++) 
		mean_abs[i] = (double *)malloc(nCols * sizeof(double));	
	abs_mean = (double **)malloc(nRows* sizeof(double *));
	for (int i = 0; i < nRows; i++) 
		abs_mean[i] = (double *)malloc(nCols * sizeof(double));	
		
	mean_freqbin_re = (double **)malloc(nRows* sizeof(double *));
	for (int i = 0; i < nRows; i++) 
		mean_freqbin_re[i] = (double *)malloc(nCols * sizeof(double));	
		
	mean_freqbin_im = (double **)malloc(nRows* sizeof(double *));
	for (int i = 0; i < nRows; i++) 
		mean_freqbin_im[i] = (double *)malloc(nCols * sizeof(double));	
		
	mean_timebin_re = (double **)malloc(nRows* sizeof(double *));
	for (int i = 0; i < nRows; i++) 
		mean_timebin_re[i] = (double *)malloc(nCols * sizeof(double));	
		
	mean_timebin_im = (double **)malloc(nRows* sizeof(double *));
	for (int i = 0; i < nRows; i++) 
		mean_timebin_im[i] = (double *)malloc(nCols * sizeof(double));	
		
	Asc_model::nRows=nRows;
    Asc_model::nCols=nCols;
	
	node::set_size(nRows,nCols);
	
}

//abs  first example
long Asc_model::abs_first()
{
	node *p;
	p= head;
	
		for (int i=0;i<nRows;i++)
			for(int j=0;j<nCols;j++)
			{
			
					first_abs[i][j]=sqrt(p->spectrogram_re[i][j]*p->spectrogram_re[i][j]+p->spectrogram_im[i][j]*p->spectrogram_im[i][j]);
			
				
			}
	
}
//mean all
long Asc_model::mean_all()
{	
    int **count;
	count = (int **)malloc(nRows* sizeof(int *));
	for (int i = 0; i < nRows; i++) 
		count[i] = (int *)malloc(nCols * sizeof(int));
	for (int i=0;i<nRows;i++)
			for(int j=0;j<nCols;j++)
			{
				count[i][j]=0;
				mean_re[i][j]=0;
				mean_im[i][j]=0;
				mean_abs[i][j]=0;
			}	
	
	node *p;
	p= head;
	while (p != NULL)
	{
		for (int i=0;i<nRows;i++)
			for(int j=0;j<nCols;j++)
			{
				if (isnan(p->spectrogram_re[i][j])==0 && isnan(p->spectrogram_im[i][j])==0)
				{
					mean_re[i][j]=mean_re[i][j]+p->spectrogram_re[i][j];
					mean_im[i][j]=mean_im[i][j]+p->spectrogram_im[i][j];
					mean_abs[i][j]=mean_abs[i][j]+sqrt(p->spectrogram_re[i][j]*p->spectrogram_re[i][j]+p->spectrogram_im[i][j]*p->spectrogram_im[i][j]);
					count[i][j]++;
				}
				
			}
		p=p->next;		
			
		
	}
    for (int i=0;i<nRows;i++)
		for(int j=0;j<nCols;j++)
		{
			mean_re[i][j]=mean_re[i][j]/count[i][j];
			mean_im[i][j]=mean_im[i][j]/count[i][j];
			mean_abs[i][j]=mean_abs[i][j]/count[i][j];
			
		}
		
	 for (int i=0;i<nRows;i++)
		for(int j=0;j<nCols;j++)
		{
			abs_mean[i][j]=sqrt(mean_re[i][j]*mean_re[i][j]+mean_im[i][j]*mean_im[i][j]);
		} 
	

    	p= head;
	for (int j=0;j<nCols;j++)
	{
		p= head;
		while (p != NULL)
		{
			
				for(int i=0;i<nRows;i++)
				{
					if (isnan(p->spectrogram_re[i][j])==0 && isnan(p->spectrogram_im[i][j])==0)
					{
						mean_timebin_re[i][j]=mean_timebin_re[i][j]+p->spectrogram_re[i][j];
						mean_timebin_im[i][j]=mean_timebin_im[i][j]+p->spectrogram_im[i][j];
						
					}
					
				}
					
			p=p->next;	
			
		}
		
		for(int i=0;j<nRows;j++)
		{
					
			mean_timebin_re[i][j]=mean_timebin_re[i][j]/count[i][j];
			mean_timebin_im[i][j]=mean_timebin_im[i][j]/count[i][j];
									
					
		}
    }	
	
	
		p= head;
	for (int i=0;i<nRows;i++)
	{
		p= head;
		while (p != NULL)
		{
			
				for(int j=0;j<nCols;j++)
				{
					if (isnan(p->spectrogram_re[i][j])==0 && isnan(p->spectrogram_im[i][j])==0)
					{
						mean_freqbin_re[i][j]=mean_freqbin_re[i][j]+p->spectrogram_re[i][j];
						mean_freqbin_im[i][j]=mean_freqbin_im[i][j]+p->spectrogram_im[i][j];
						
					}
					
				}
					
			p=p->next;	
			
		}
		
		for(int j=0;j<nRows;j++)
		{
					
			mean_freqbin_re[i][j]=mean_freqbin_re[i][j]/count[i][j];
			mean_freqbin_im[i][j]=mean_freqbin_im[i][j]/count[i][j];
									
					
		}
    }
	
	for (int i = 0; i < nRows; i++) 
		free(count[i]);
    free(count);	

	 
}

long Asc_model::save(char *mname_a)
{
    FILE* fp = fopen(mname_a, "w");
	 if (fp == NULL) {
    fprintf(stdout, "*> Asc_model::save: error creating file (%s)\n",
	    mname_a);
    return Asc_model::ERR_FILE;
    }
	node *p;
	p= head;
	fwrite(&nRows, sizeof(long), 1, fp);
	fwrite(&nCols, sizeof(long), 1, fp);
	
	fwrite(&node::len, sizeof(long), 1, fp);
	
	fwrite(&A1_score1,sizeof(double),1,fp);
	fwrite(&B1_score1,sizeof(double),1,fp);
	fwrite(&A1_score2,sizeof(double),1,fp);
	fwrite(&B1_score2,sizeof(double),1,fp);
	fwrite(&A1_score3,sizeof(double),1,fp);
	fwrite(&B1_score3,sizeof(double),1,fp);
	
	while (p != NULL)
	{
	  
				
	  fwrite(p->filename,sizeof(char),1000,fp);
      fwrite(&p->example_no,sizeof(long),1,fp);	  
	  for (int i = 0; i < nRows; i++) 
	  fwrite(p->spectrogram_re[i], sizeof(double),nCols,fp);
	  for (int i = 0; i < nRows; i++) 
	  fwrite(p->spectrogram_im[i], sizeof(double),nCols,fp); 
		p=p->next;
	}
	
	fclose(fp);
	return Asc_model::NO_ERR;
}

long Asc_model::load(char *mname_a)
{
    FILE* fp = fopen(mname_a, "r");
	double **sp_re=NULL;
	double **sp_im=NULL;
	char filename[1000];
	long example_no;
	
	long len;
	long count=0;
	 if (fp == NULL) {
    fprintf(stdout, "*> Asc_model::load: error reading file (%s)\n",
	    mname_a);
    return Asc_model::ERR_FILE;
    }
	
	fread(&nRows, sizeof(long), 1, fp);
	fread(&nCols, sizeof(long), 1, fp);
	init(nRows,nCols);
	
	fread(&len, sizeof(long), 1, fp);
	
	fread(&A1_score1,sizeof(double),1,fp);
	fread(&B1_score1,sizeof(double),1,fp);
	fread(&A1_score2,sizeof(double),1,fp);
	fread(&B1_score2,sizeof(double),1,fp);
	fread(&A1_score3,sizeof(double),1,fp);
	fread(&B1_score3,sizeof(double),1,fp);
	 sp_re = (double **)malloc(nRows * sizeof(double *));
	for (int i = 0; i < nRows; i++) 
		sp_re[i] = (double *)malloc(nCols * sizeof(double));
	sp_im = (double **)malloc(nRows* sizeof(double *));
	for (int i = 0; i < nRows; i++) 
		sp_im[i] = (double *)malloc(nCols * sizeof(double));

		while (count<len)
		{
		  
		  fread(filename,sizeof(char),1000,fp);
		  fread(&example_no,sizeof(long),1,fp);	  
		 
		  for(int i=0;i<nRows;i++)
		  fread(sp_re[i], sizeof(double),nCols,fp);
		   for(int i=0;i<nRows;i++)
		  fread(sp_im[i], sizeof(double),nCols,fp);
		  
		  add_example(sp_re,sp_im,nRows,nCols,filename);
		   
			count++;
		}
		
	
	
	fclose(fp);

	for (int i = 0; i < nRows; i++) 
		free(sp_re[i]);
    free(sp_re);
	
	for (int i = 0; i < nRows; i++) 
		free(sp_im[i]);
    free(sp_im);
	
	return Asc_model::NO_ERR;
}

long Asc_model::print()
{
	node *p;
	p=head;
	while(p!=NULL)
	{
		printf("\n%ld :%s\n",p->example_no,p->filename);
		printf("spectrogram_re\n");
		for (int i=0;i<nRows;i++)
		{
			for(int j=0;j<nCols;j++)
				printf("%f ",p->spectrogram_re[i][j]);
			printf("\n");
		}
		
		printf("spectrogram_im\n");
		for (int i=0;i<nRows;i++)
		{
			for(int j=0;j<nCols;j++)
				printf("%f ",p->spectrogram_im[i][j]);
			printf("\n");
		}
		
		
		p=p->next;	
	}	

    return 0;
}
long Asc_model::retrieve()
{
  node *p;
  p=head;
  long count=1;
  while (p!=NULL)
  {
	fprintf(stdout,"record %ld :\n%s\n",count,p->filename);
	count++;
	p=p->next;
  }
  
  return 0;
}