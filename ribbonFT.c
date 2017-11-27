#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include <stdlib.h>
#include<strings.h>
#include<stdarg.h>//Needed for the function print_and_exit()
#include "fftw3.h"
#define M_PI 3.14159265358979323846

#define MAXLENGTH 10000000
#define PERIOD 10000
#define RUNMAX 201
#define NMAX 2000
#define MAXFRAMES 20001

double v[MAXLENGTH]; // The list of measures
double h_width[MAXFRAMES][NMAX];
double power_spectrum_frame[MAXFRAMES][NMAX];
double avg_hFT[RUNMAX][NMAX];
double avg_h[RUNMAX][NMAX];
double hFT_width_avg[NMAX];//Averaging fourier amplitude square across runs
double jk_blocks[NMAX][RUNMAX];
double jk_error[NMAX];
double error[NMAX]; //RMSE error in |hFT|^2

int NX,NY,RUNS,STEPS,FRAMES,JK_BIN_COUNT;
double KAPPA;
double NORM;

void print_and_exit(char *, ...); //Print out an error message and exits

 
int main(int argc, char **argv)
{

   FILE *Fin,*fp;
   char data_file[1024];
   int i,j,n;
   double *hd;
   fftw_complex *hFT;
   fftw_plan pdir;

   switch (argc){
   case 5:
       NX = atoi(argv[1]);
       NY = atoi(argv[2]);
       KAPPA = atof(argv[3]);
       STEPS = atoi(argv[4]);
       break;
   default:
       print_and_exit("Usage Pass command line arguments:NX NY KAPPA STEPS\n");
   }
 
   FRAMES = STEPS/PERIOD;
   //printf("frames %d\n",frames);
   n = NX; // Number of real data points for which FFT is taken

   char validrunfile[1024],validrunpath[1024];
   int nx,total_frames;//reading from file
   RUNS = 0;
  
   sprintf(validrunfile,"../Sim_dump_ribbon/L%d/W%d/k%.1f/valid_thermal_runs.log",NX,NY,KAPPA);

   if(NULL==(fp=fopen(validrunfile,"r")))
        print_and_exit("I could not open file with simulation valid run numbers %s\n",validrunfile);

   /*  Allocating memory for FFT       */
   hd = (double *) fftw_malloc(sizeof(double)*n);
   hFT = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*((n/2)+1));

   // Plan for FFTW    
   pdir = fftw_plan_dft_r2c_1d(n,hd,hFT,FFTW_MEASURE);
 
   while (fscanf(fp, "%s", validrunpath) == 1)// 1 is returned if fscanf reads a valid run path
   {
	   //fscanf (file, "%d", &runnum); // file contains the runs to be analyzed
	   sprintf(data_file,"%s/hgt_widthavg.bin",validrunpath);
	   //printf("%s\n",data_file);
	   if(NULL==(Fin=fopen(data_file,"rb")))
	      print_and_exit("I could not open binary file %s\n",data_file);

	    //We read the data file
	    //fread(h_width,sizeof(double),MAXFRAMES*NMAX,Fin);//frames*NX,Fin);
	    fread(&nx,sizeof(int),1,Fin);
  	    fread(&total_frames,sizeof(int),1,Fin);
	    for(int i=0;i<total_frames;i++)
  	    {
        	for(int j=0;j<nx;j++)
        	{
                	fread(&h_width[i][j],sizeof(double),1,Fin);
        	}
  	    }
	    fclose(Fin); 

	    if(nx != NX || total_frames != FRAMES)
		print_and_exit("Mismatch between nx or total frames read from bin file and calculated values\n");

	    // Now fill in the vector hd; initializing of input array should be done after creating the plan
	    for(j=FRAMES/2;j<=FRAMES;j++)
	    {
		//printf("%d\n",j);
		for(i=0; i<n; i++)
		{
			hd[i] = h_width[j][i];
		}

	   	//Execute the FFTW
	   	fftw_execute(pdir);

	   	// hFT contains the FT of hd, we need to compute | hFT | ^ 2
	   	for (i=0;i<((n/2)+1);i++)
	   	{
			//hFT[i][0] = (hFT[i][0]*hFT[i][0] + hFT[i][1]*hFT[i][1])/(pow(signal_length/2,2));
			//hFT[i][1] = 0;
			power_spectrum_frame[j-(FRAMES/2)][i] = (hFT[i][0]*hFT[i][0] + hFT[i][1]*hFT[i][1])/(pow(n/2,2));
	   	}
	     }

	     //Adding fourier amplitudes at same x for different frames same run
     for(i=0;i<((n/2)+1);i++)
	{
		avg_hFT[RUNS][i]=0;
		for(j=0;j<FRAMES/2;j++)
		{
			avg_hFT[RUNS][i]+=power_spectrum_frame[j][i];
		}
		avg_hFT[RUNS][i]/=(FRAMES/2);
		NORM = 1;//avg_hFT[RUNS][1];
		//printf("%d\t%.8f\n",run,avg_hFT[run][i]);
	}

	for(i=0;i<((n/2)+1);i++)
        {
                avg_hFT[RUNS][i]/=NORM;
                //printf("%d\t%.8f\n",RUNS,avg_hFT[RUNS][i]);
        }	

	RUNS++; 
    }

    fclose(fp);

    JK_BIN_COUNT = RUNS; //Number of Jack Knife bins
	   

    /* Average of power spectrum across runs	*/
    for(i=0;i<((n/2)+1);i++)
    {
	hFT_width_avg[i]=0;
	//printf("%d",i);
	for(int r=0;r<RUNS;r++)
	{
		hFT_width_avg[i]+=avg_hFT[r][i];
	}
	hFT_width_avg[i]/=RUNS;
	//printf("%.8f",hFT_width_avg[i]);
    }

   /*      RMSE error in power spectrum   */
    for(i=0;i<((n/2)+1);i++)
    {
	//printf("%d\t%.8g\t",i,hFT_width_avg[i]);
	//if (i==1)
		//printf("%.8f\t%.8f\n",sqrt(3)*(NY-1)/(2*NX-1),hFT_width_avg[i]/hFT_width_avg[i+1]);
	error[i]=0;
	for(int r=0;r<RUNS;r++)
	{
		error[i] += pow((avg_hFT[r][i] - hFT_width_avg[i]),2);
	}
	error[i] = sqrt(error[i]/RUNS);
	//printf("%.8g\n",error[i]);
    }  

    /*		Jack Knife Error estimation	*/
    
    /*		Total of Jack Knife blocks at each sampling interval	*/
    for(i=0;i<((n/2)+1);i++)
    {
	jk_blocks[i][JK_BIN_COUNT]=0;
	for(j=0;j<JK_BIN_COUNT;j++)
	{
		jk_blocks[i][JK_BIN_COUNT] += avg_hFT[j][i]; //summing avg_hFT for all runs at each i
	}
    }    

    /*		Jack Knife Blocking	*/
    for(i=0;i<((n/2)+1);i++)
    {
        for(j=0;j<JK_BIN_COUNT;j++)
        { 
		jk_blocks[i][j]= (jk_blocks[i][JK_BIN_COUNT] - avg_hFT[j][i])/(JK_BIN_COUNT-1);	
	}
    }

    double jk_error_term1[NMAX],jk_error_term2[NMAX];
    /*		Jack Knife Error	*/
    for(i=0;i<((n/2)+1);i++)
    {
	jk_error[i]=0;
	jk_error_term1[i]=0;
	jk_error_term2[i]=0;
    } 

    for(i=0;i<((n/2)+1);i++)
    {
        for(j=0;j<JK_BIN_COUNT;j++)
        {
		jk_error_term1[i] += jk_blocks[i][j] * jk_blocks[i][j];
	}
	jk_error_term1[i] = (1.0/JK_BIN_COUNT) * jk_error_term1[i];
    }

    for(i=0;i<((n/2)+1);i++)
    {
        for(j=0;j<JK_BIN_COUNT;j++)
        {
		jk_error_term2[i] += (1.0/JK_BIN_COUNT) * jk_blocks[i][j];
	}
	jk_error_term2[i] = jk_error_term2[i] * jk_error_term2[i];
        /*	JK Error	*/
	jk_error[i] = sqrt((JK_BIN_COUNT-1)*(jk_error_term1[i] - jk_error_term2[i]));
	printf ("%d\t%.8g\t%.8g\t%.8g\n",i,i*(2*M_PI/n),jk_blocks[i][JK_BIN_COUNT]/JK_BIN_COUNT,jk_error[i]);	
    }

    //printf("Cleaning up\n");
    fftw_destroy_plan(pdir);
    fftw_free(hd);
    fftw_free(hFT);
    return 0;   
}

void print_and_exit(char *format, ...)
{
    va_list list;

    va_start(list,format);
    vprintf(format,list);
    va_end(list);
    exit(1);
}

