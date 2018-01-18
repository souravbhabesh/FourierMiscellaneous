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
int filter_survival_time[MAXLENGTH],filter_orientation[MAXLENGTH];

int NX,NY,RUNS,STEPS,FRAMES;
double KAPPA;

void print_and_exit(char *, ...); //Print out an error message and exits
int left_shift_array(int [], int , int );
 
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

   //Initializing the output file
   sprintf(out_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/survivaltime.log",NX,NY,KAPPA);
   printf("Output file: %s\n",out_file);
   if(NULL==(Fout=fopen(out_file,"w")))
         print_and_exit("I could not open file with simulation survival time %s\n",out_file);

   //We read the data file
   fread(&RUNNUM,sizeof(int),1,Fin);
   printf("Total Runs: %d\n",RUNNUM);
   RUNNUM=1;//reading just the 2nd run for checking 
   for(int i=0;i<RUNNUM;i++)
   {
        for(int j=0;j<FRAMES;j++)
        {
                fread(&cnode[i][j],sizeof(double),1,Fin);
        }
   }
   fclose(Fin);

   int cnt,time_up,time_down,time_tran,threshold=2,scnt=0; 
   for(i=0;i<RUNNUM;i++)
   {
        time_up=0;
        time_down=0;
        time_tran=0;
        cnt=0;
        //fprintf(Fout,"%d\t%.8f\t%d\t%d\t%d\n",998,cnode[i][998],orientation[cnt],threshold,-1*threshold); 

	for(int j=0.2*FRAMES;j<=FRAMES;j++)//Discard first 20% data
        {
                // orientation: Up (+1) Transition (0) Down (-1)
		if(cnode[i][j]<-2)  //Down
		{
                  if(cnode[i][j-1]>-2)//cross over: transition to Down
                  {
                    survival_time[cnt]=time_tran;
                    //orientation[cnt]=0;//  -1:Down 0:Transition 1:Up
                    time_tran=0;
                    cnt++;
                  }
                  time_down++;
                  orientation[cnt]=-1;
                }

                if (cnode[i][j]>2)  //Up
                {
                  if(cnode[i][j-1]<2) //cross over: transition to Up
                  {
                    survival_time[cnt]=time_tran;
                    time_tran=0;
                    cnt++;
                  }
                  time_up++;
                  orientation[cnt]=1;
		}

                if(cnode[i][j]>-2 && cnode[i][j]<2) // In Transition
                {
                  if(cnode[i][j-1]<-2) // cross over: Down to Transition
                  {
                    survival_time[cnt]=time_down;
                    time_down=0;
                    cnt++;
                  }

                  if(cnode[i][j-1]>2) // cross over : Up to Transition
                  {
                    survival_time[cnt]=time_up;
                    time_up=0;
                    cnt++;
                  }
                  time_tran++;
                  orientation[cnt]=0;
                }
            //Printing to output file
            fprintf(Fout,"%d\t%.8f\t%d\t%d\t%d\t%d\t%d\n",j,cnode[i][j],orientation[cnt],threshold,-1*threshold,cnt,survival_time[cnt]);
        }
    

    // Eliminating small excursions to Transition from Up/Down
    for (int k=1;k<cnt;k++)
    {
      if(((orientation[k-1]==-1 && orientation[k+1]==-1 && survival_time[k-1]>10) || (orientation[k-1]==1 && orientation[k+1]==1 && survival_time[k-1]>10)) && orientation[k]==0)
      {
          //printf("Excursion Survival time %d : %d\n",i,survival_time[k]);
          if(survival_time[k]<10)
          {
              survival_time[k-1]+=survival_time[k]+survival_time[k+1];
              left_shift_array(survival_time,cnt-2,k);
              left_shift_array(orientation,cnt-2,k);
              cnt-=2;//since we have shifted two places left
              k--;//Next k value remains same due to shift
          }
      }
    }
/*
    // Eliminating small excursions to Up/Down from Transition
    for (int k=1;k<cnt;k++)
    {
        if(orientation[k-1]==0 && orientation[k+1]==0 && (orientation[k]==1  || orientation[k]==-1))
        {
            //printf("Survival time %d : %d at %d orientation %d\n",k,survival_time[k],j,orientation[k]);
            if(survival_time[k]<10)
            {
              survival_time[k-1]+=survival_time[k]+survival_time[k+1];
              left_shift_array(survival_time,cnt-2,k);
              left_shift_array(orientation,cnt-2,k);
              cnt-=2;//since we have shifted two places left
              k--;//Next k value remains same due to shift
            }
        }
    }


    // Eliminating small excursions to Transition from Up/Down
    for (int k=1;k<cnt;k++)
    {
      if(((orientation[k-1]==-1 && orientation[k+1]==-1 && survival_time[k-1]>10) || (orientation[k-1]==1 && orientation[k+1]==1 && survival_time[k-1]>10)) && orientation[k]==0)
      {
          //printf("Excursion Survival time %d : %d\n",i,survival_time[k]);
          if(survival_time[k]<10)
          {
              survival_time[k-1]+=survival_time[k]+survival_time[k+1];
              left_shift_array(survival_time,cnt-2,k);
              left_shift_array(orientation,cnt-2,k);
              cnt-=2;//since we have shifted two places left
              k--;//Next k value remains same due to shift
          }
      }
    }

    // Eliminating small excursions to Up/Down from Transition
    for (int k=1;k<cnt;k++)
    {
        if(orientation[k-1]==0 && orientation[k+1]==0 && (orientation[k]==1  || orientation[k]==-1))
        {
            //printf("Survival time %d : %d at %d orientation %d\n",k,survival_time[k],j,orientation[k]);
            if(survival_time[k]<10)
            {
              survival_time[k-1]+=survival_time[k]+survival_time[k+1];
              left_shift_array(survival_time,cnt-2,k);
              left_shift_array(orientation,cnt-2,k);
              cnt-=2;//since we have shifted two places left
              k--;//Next k value remains same due to shift
            }
        }
    }
*/
    //printf("cnt %d\n",cnt);
    int time_cnt=0.2*FRAMES;
    for(int k=0;k<=cnt;k++)
    {
        /*
        if(orientation[k]!=0)
        {
            printf("%d\t%d\t%d\n",scnt,survival_time[k],orientation[k]);
            //printf("%d\n",survival_time[k]);
            scnt++;
        }
        */
        
        for(int l=0; l < survival_time[k];l++)
        {
            printf("%d\t%d\t%d\t%d\n",time_cnt,orientation[k],survival_time[k],k);
            time_cnt++;
            
        }
        
    }

   }

fclose(Fout);
return 0;   
}

//Left shift array by 2 places
int left_shift_array(int arr[], int size, int site)
{
    for(int i=site;i<size;i++)
    {
        arr[i]=arr[i+2];
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

