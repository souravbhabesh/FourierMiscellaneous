#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include <stdlib.h>
#include<strings.h>
#include<stdarg.h>//Needed for the function print_and_exit()
//#define M_PI 3.14159265358979323846

#define MAXLENGTH 10000000
#define PERIOD 10000
#define RUNMAX 201
#define NMAX 2000
#define MAXFRAMES 20001

double cnode[RUNMAX][MAXFRAMES];

int NX,NY,RUNS,STEPS,FRAMES;
double KAPPA;

void print_and_exit(char *, ...); //Print out an error message and exits

 
int main(int argc, char **argv)
{

   FILE *Fin;
   char data_file[1024];
   int i,RUNNUM,tot_runs;

   switch (argc){
   case 6:
       NX = atoi(argv[1]);
       NY = atoi(argv[2]);
       KAPPA = atof(argv[3]);
       STEPS = atoi(argv[4]);
       RUNNUM = atoi(argv[5]);
       break;
   default:
       print_and_exit("Usage Pass command line arguments:NX NY Kappa steps RUNNUM\n");
   }

   FRAMES = STEPS/PERIOD;

   sprintf(data_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/cnode.bin",NX,NY,KAPPA);
   //printf("%s\n",data_file);
   if(NULL==(Fin=fopen(data_file,"rb")))
   	print_and_exit("I could not open binary file %s\n",data_file);

   //printf("frames %d\n",frames);

   //We read the data file
   //fread(cnode,sizeof(double),RUNMAX*MAXFRAMES,Fin);
   fread(&tot_runs,sizeof(int),1,Fin);//Number of total runs
   for(int i=0;i<tot_runs;i++)
   {
        for(int j=0;j<FRAMES;j++)
        {
                fread(&cnode[i][j],sizeof(double),1,Fin);
        }
   }
   fclose(Fin);

   if(RUNNUM > tot_runs)
	print_and_exit("RUNNUM input is greater than total runs \n");

   for(i=0;i<FRAMES;i++)
   {
	printf("%d\t%.8f\n",i*PERIOD,cnode[RUNNUM-1][i]);
   }
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

