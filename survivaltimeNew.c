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

double cnode[RUNMAX][MAXFRAMES];
int survival_time[MAXLENGTH],orientation[MAXLENGTH];

int NX,NY,RUNS,STEPS,FRAMES;
double KAPPA;

void print_and_exit(char *, ...); //Print out an error message and exits

 
int main(int argc, char **argv)
{

   FILE *Fin,*Fout;
   char data_file[1024],out_file[256];
   int i,RUNNUM;

   switch (argc){
   case 5:
       NX = atoi(argv[1]);
       NY = atoi(argv[2]);
       KAPPA = atof(argv[3]);
       STEPS = atoi(argv[4]);
       break;
   default:
       print_and_exit("Usage Pass command line arguments:NX NY Kappa STEPS\n");
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

   int start_time = -1, count=0,time_up,time_down; 
   for(i=0;i<RUNNUM;i++)
   {
        time_up=0;
        time_down=0;
        time_tran=0;
	for(int j=0.2*FRAMES;j<FRAMES;j++)//Discard first 20% data
        {
		if(cnode[i][j]<-2)
		{
		  time_down++;
                  if(cnode[i][j-1]>-2)
                  {
                    survival_time[cnt]=time_tran;
                    orientation[cnt]=2;//0-Down 1-Up 2- in transition
                    time_tran=0;
                    cnt++;
                  }

                }
                if (cnode[i][j]>2)
                {
                  time_up++; 
                  if(cnode[i][j-1]<2)
                  {
                    survival_time[cnt]=time_tran;
                    orientation[cnt]=2;//0-Down 1-Up 2- in transition
                    time_tran=0;
                    cnt++;
                  }

		}
                if(cnode[i][j]>-2 || cnode[i][j]<2)
                {
                  time_tran++;
                  if(cnode[i][j-1]<-2)
                  {
                    survival_time[cnt]=time_down;
                    orientation[cnt]=0;//0-Down 1-Up
                    time_down=0;
                    cnt++;
                  }

                  if(cnode[i][j-1]>2)
                  {
                    survival_time[cnt]=time_up;
                    orientation[cnt]=1;//0-Down 1-Up
                    time_up=0;
                    cnt++;
                  }
                }
	}
	
   }


  sprintf(out_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/survivaltime.log",NX,NY,KAPPA);
  printf("Output file: %s\n",out_file);
  if(NULL==(Fout=fopen(out_file,"w")))
        print_and_exit("I could not open file with simulation survival time %s\n",out_file);
  for (int i=0;i<count;i++)
  {
        fprintf(Fout,"%d\n",survival_time[i]);
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

