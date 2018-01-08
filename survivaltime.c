#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include <stdlib.h>
#include<strings.h>
#include<stdarg.h>//Needed for the function print_and_exit()
#define M_PI 3.14159265358979323846

#define MAXLENGTH 10000000
#define PERIOD 10000
#define RUNMAX 201
#define NMAX 2000
#define MAXFRAMES 20001
#define MAXBINS 100

double cnode[RUNMAX][MAXFRAMES];
int survival_time[MAXLENGTH];

int NX,NY,RUNS,STEPS,FRAMES,CUTOFF;
double KAPPA;

void print_and_exit(char *, ...); //Print out an error message and exits

 
int main(int argc, char **argv)
{

   FILE *Fin,*Fout;
   char data_file[1024],out_file[256];
   int i,RUNNUM;

   switch (argc){
   case 6:
       NX = atoi(argv[1]);
       NY = atoi(argv[2]);
       KAPPA = atof(argv[3]);
       STEPS = atoi(argv[4]);
       CUTOFF = atoi(argv[5]);
       break;
   default:
       print_and_exit("Usage Pass command line arguments:NX NY Kappa STEPS (SURVIVAL TIME CUTOFF)\n");
   }

   FRAMES = STEPS/PERIOD;

   sprintf(data_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/cnode.bin",NX,NY,KAPPA);
   printf("%s\n",data_file);
   if(NULL==(Fin=fopen(data_file,"rb")))
   	print_and_exit("I could not open binary file %s\n",data_file);
   

   //printf("frames %d\n",frames);

   //We read the data file
   fread(&RUNNUM,sizeof(int),1,Fin);
   for(int i=0;i<RUNNUM;i++)
   {
        for(int j=0;j<FRAMES;j++)
        {
                fread(&cnode[i][j],sizeof(double),1,Fin);
        }
   }
   fclose(Fin);

   int start_time = -1, count=0,time; 
   for(i=0;i<RUNNUM;i++)
   {
	for(int j=0.2*FRAMES;j<FRAMES;j++)//Discard first 20% data
        {
		if((cnode[i][j-1]<0 && cnode[i][j]>0) || (cnode[i][j-1]>0 && cnode[i][j]<0))
		{
			time = j-1-start_time;
			if(start_time != -1 && time != 0 && time > CUTOFF)//time=0 cross over back to back
                        {
                                survival_time[count]=time;
				//if(time == 0)
					//printf("%d %d\n",j-1,start_time);
                                count++;
                        }
                        start_time = j;
		}
	}
	start_time = -1;
   }


  sprintf(out_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/survivaltime.log",NX,NY,KAPPA);
  printf("Output file: %s\n",out_file);
  if(NULL==(Fout=fopen(out_file,"w")))
        print_and_exit("I could not open file with simulation survival time %s\n",out_file);

/*
  for (int i=0;i<count;i++)
  {
        fprintf(Fout,"%d\n",survival_time[i]);
  }
*/

  int freq[MAXBINS];
  //Binning data to make histogram
  for(int i=1;i<MAXBINS;i++)
  {
    freq[i-1]=0;
    for(int j=0;j<count;j++)
    {
      if((survival_time[j] > (i*CUTOFF)) &&  ((survival_time[j] < (i+1)*CUTOFF)))
        freq[i-1]++;
    }
    if(freq[i-1]>0)
      fprintf(Fout,"%d\t%d\n",i,freq[i-1]);
  }
  fclose(Fout);
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

