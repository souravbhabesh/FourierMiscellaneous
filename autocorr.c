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
#define NMAX 20000
#define MAXFRAMES 20001

double cnode[RUNMAX][MAXFRAMES];
double data[RUNMAX][MAXFRAMES];
double autocorr[RUNMAX][MAXFRAMES];
double norm[RUNMAX];
double jk_error[NMAX];
double error[NMAX]; //RMSE error in |hFT|^2

int NX,NY,RUNS,STEPS,FRAMES,JK_BIN_COUNT;
double KAPPA;

void print_and_exit(char *, ...); //Print out an error message and exits

 
int main(int argc, char **argv)
{

   FILE *Fin;
   char data_file[1024];
   int i,j,n;
   double *hd,*corr_h;
   fftw_complex *hFT;
   fftw_plan pdir,pinv;

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
   n = 0.8*FRAMES; // Number of real data points for which FFT is taken, discarding top 20% data
   //printf("n %d\n",n);
 
   sprintf(data_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/cnode.bin",NX,NY,KAPPA);
   //printf("%s\n",data_file);
   if(NULL==(Fin=fopen(data_file,"rb")))
   	print_and_exit("I could not open binary file %s\n",data_file);


   //We read the data file
   fread(&RUNS,sizeof(int),1,Fin);
   //RUNS=10;
   //printf("RUNS %d\n",RUNS);
   for(int i=0;i<RUNS;i++)
   {
        for(int j=0;j<FRAMES;j++)
        {
                fread(&cnode[i][j],sizeof(double),1,Fin);
        }
   }
   fclose(Fin);

   // Allocating space for carrying out FFT's
   hd = (double *) fftw_malloc(sizeof(double)*n*2);
   corr_h = (double *) fftw_malloc(sizeof(double)*n*2); 
   hFT = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*(n+1));
             
    // We prepare the FFTW
    pdir = fftw_plan_dft_r2c_1d(2*n,hd,hFT,FFTW_MEASURE);
    pinv = fftw_plan_dft_c2r_1d(2*n,hFT,corr_h,FFTW_MEASURE);

   double avg;
   int iblo,start=0.2*FRAMES;
   //printf("start %d\n",start);
   for(i=0; i<n; i++)
   {
     for(iblo=0; iblo<RUNS; iblo++) // FIX NUMBER OF RUNS
       data[iblo][i] = cnode[iblo][start+i];
   }
     	
   
   for(iblo=0; iblo<RUNS; iblo++)
   {
     	   avg=0;//for each run

	   for(i=0; i<n; i++)
	     avg += data[iblo][i];

	   avg /= n;

	   // Now fill in the vector hd
	   for(i=0; i<n; i++)
	     hd[i] = data[iblo][i] - avg;
	   for(i=n; i<2*n; i++)
	     hd[i]=0;

	   // We execute the FFTW
	   fftw_execute(pdir); 

	   // hFT contains the FT of h, we need to compute | hFT | ^ 2
           for (i=0;i<n+1;i++){
      		hFT[i][0] = hFT[i][0]*hFT[i][0] + hFT[i][1]*hFT[i][1];
      		hFT[i][1] = 0;
    		}
    	   fftw_execute(pinv);

	   // The inverse FT of | hFT | ^ 2 is autocorrelation, but we must
    	   //normalize
    	   for(i=0; i<n; i++)
	   {
      		corr_h[i] /= (n-i);
		autocorr[iblo][i] = corr_h[i];
	   }

	   //for (i=0;i<n;i++)
      		//printf("%d %.8g\n", i, (hFT[i][0])/(pow(n,2)));

	  //for (i=0;i<n;i++)
      		//printf("%d %.8g\n", i, corr_h[i]/corr_h[0]); 
   }

/*
   for (i=0;i<n;i++)
   	printf("%d %.8g %.8g %.8g %.8g %.8g %.8g\n", i, autocorr[0][i],autocorr[1][i],autocorr[2][i],autocorr[RUNS-2][i],autocorr[RUNS-1][i],autocorr[RUNS][i]);
*/

/*
   int f=100;
   for(i=0;i<=RUNS;i++)
	printf("%d %.8f\n",i,autocorr[i][f]);
*/

   for(i=0; i<n; i++)
   {
     autocorr[RUNS][i]=0;
     for(iblo=0; iblo<RUNS; iblo++)
	     autocorr[RUNS][i]+= autocorr[iblo][i]; 
     for(iblo=0; iblo<RUNS; iblo++)
	     autocorr[iblo][i] = (autocorr[RUNS][i]- autocorr[iblo][i])/(RUNS-1);
     autocorr[RUNS][i] /= RUNS;
   }
     for(iblo=0; iblo<=RUNS; iblo++)
	     norm[iblo] = autocorr[iblo][0];
   for(i=0; i<n; i++)
   {
     for(iblo=0; iblo<=RUNS; iblo++)
	     autocorr[iblo][i]/= norm[iblo];
   }



    double jk_error_term1[NMAX],jk_error_term2[NMAX];
    //		Jack Knife Error	
    for(i=0;i<n;i++)
    {
	jk_error[i]=0;
	jk_error_term1[i]=0;
	jk_error_term2[i]=0;
    } 

    JK_BIN_COUNT=RUNS; //FIX RUN NUMBERS
    for(i=0;i<n;i++)
    {
        for(j=0;j<JK_BIN_COUNT;j++)
        {
		jk_error_term1[i] += (autocorr[j][i] * autocorr[j][i]);
	}
	jk_error_term1[i] = (1.0/JK_BIN_COUNT) * jk_error_term1[i];
    }

    for(i=0;i<n;i++)
    {
        for(j=0;j<JK_BIN_COUNT;j++)
        {
		jk_error_term2[i] += (1.0/JK_BIN_COUNT) * autocorr[j][i];
	}
	jk_error_term2[i] = jk_error_term2[i] * jk_error_term2[i];
        //	JK Error	
	jk_error[i] = sqrt(((JK_BIN_COUNT-1)*(jk_error_term1[i]-jk_error_term2[i])));
	printf ("%d\t%.8g\t%.8g\n",i,autocorr[RUNS][i],jk_error[i]);
	//printf ("%d\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\n",i,autocorr[JK_BIN_COUNT][i],jk_error[i],jk_error_term1[i],jk_error_term2[i],sqrt(jk_error_term1[i]-jk_error_term2[i]));	
    }

   fftw_destroy_plan(pdir);
   fftw_destroy_plan(pinv);
   fftw_free(hd);
   fftw_free(hFT);
   fftw_free(corr_h);

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

