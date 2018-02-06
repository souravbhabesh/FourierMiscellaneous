#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include <stdlib.h>
#include<strings.h>
#include<stdarg.h>//Needed for the function print_and_exit()
#include "nrutil.h"
#include "nrutil.c"

#define M_PI 3.14159265358979323846
#define SWAP(a,b); {temp=(a);(a)=(b);(b)=temp;}
#define MAXLENGTH 10000
#define PERIOD 10000
#define RUNMAX 201
#define NMAX 2000
#define MAXFRAMES 20001
#define M 7
#define NSTACK 60

double cnode[RUNMAX][MAXFRAMES];
int survival_time[MAXLENGTH],orientation[MAXLENGTH],run_id[MAXLENGTH];
int filter_survival_time[MAXLENGTH],filter_orientation[MAXLENGTH];
int survival_time_sorted[MAXLENGTH],survival_time_unsorted[MAXLENGTH];

int NX,NY,RUNS,STEPS,FRAMES,BIN_SIZE;
double KAPPA,FILTER1,FILTER2;

void print_and_exit(char *, ...); //Print out an error message and exits
int left_shift_array(int [], int , int );
void sort(int , int []);

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
       BIN_SIZE = atoi(argv[5]);
       break;
   default:
       print_and_exit("Usage Pass command line arguments:NX NY Kappa STEPS BIN_SIZE\n");
   }

   FRAMES = STEPS/PERIOD;

   sprintf(data_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/cnode.bin",NX,NY,KAPPA);
   printf("%s\n",data_file);
   if(NULL==(Fin=fopen(data_file,"rb")))
   	print_and_exit("I could not open binary file %s\n",data_file);
   

   //printf("frames %d\n",frames);

   //Initializing the output file
   sprintf(out_file,"../Sim_dump_stretched/L%d/W%d/k%.1f/survivaltimeNew.log",NX,NY,KAPPA);
   printf("Output file: %s\n",out_file);
   if(NULL==(Fout=fopen(out_file,"w")))
         print_and_exit("I could not open file with simulation survival time %s\n",out_file);

   //We read the data file
   fread(&RUNNUM,sizeof(int),1,Fin);
   printf("Total Runs: %d\n",RUNNUM);
   //RUNNUM=1;//reading just the 1st run for checking 
   for(int i=0;i<RUNNUM;i++)
   {
        for(int j=0;j<FRAMES;j++)
        {
                fread(&cnode[i][j],sizeof(double),1,Fin);
        }
   }
   fclose(Fin);

   int j,cnt=0,time,time_down,threshold=2,scnt=0; 
   for(i=0;i<RUNNUM;i++)
   {
        j=0;
       int state=0;
        while(fabs(cnode[i][j])<2 || j<0.1*FRAMES){
                j++;
                fprintf(Fout,"%d\t%d\t%.8f\t%d\n",i,j,cnode[i][j],state);
        }

        state=(cnode[i][j]>2)?1:-1;
        orientation[cnt]=state;

        time=0;
	for(j=j+1;j<=FRAMES;j++)
        {
		if((state==1 && cnode[i][j]>-2) || (state==-1 && cnode[i][j]<2))  //Up
		{
                    time++;
                }
                else
                {
                    survival_time[cnt]=time;
                    run_id[cnt]=i;
                    time=0;
                    cnt++;
                    state*=-1;
                    orientation[cnt]=state;
                }
                fprintf(Fout,"%d\t%d\t%.8f\t%d\n",i,j,cnode[i][j],state);
          }
        }

/*
   //Printing orientation of the time series
    int time_cnt=0.2*FRAMES;
    for(int k=0;k<=cnt;k++)
    {

       for(int l=0; l < survival_time[k];l++)
        {
            printf("%d\t%d\t%d\t%d\n",time_cnt,orientation[k],survival_time[k],k);
            time_cnt++;
        }
    }
*/
/*
    printf("cnt %d\n",cnt);
    time_cnt=0.2*FRAMES;
    for(int k=0;k<=cnt;k++)
    {
            survival_time_unsorted[k]=survival_time[k]; //Unsorted survival times
            survival_time_sorted[k]=survival_time_unsorted[scnt];//This array will be sorted
            //printf("%d\t%d\t%d\n",scnt,survival_time_unsorted[scnt],survival_time_sorted[scnt]);
    }


    //Sort the survival times
    sort(scnt,survival_time_sorted);
    printf("scnt = %d\n",scnt);    
*/
    for(int k=0;k<=cnt;k++)
    {
        //Printing to check sorting is working correctly
        //printf("%d\t%d\t%d\n",k,survival_time[k],run_id[k]);
    }

    //Bin the sorted survival times
    //int BIN_SIZE=5;
    double bin[MAXLENGTH],cumulative[MAXLENGTH];
    for(int k=0;k<=MAXLENGTH;k++)
        bin[k]=0;
    for(int k=0;k<=cnt;k++)
        bin[(int)(survival_time[k]/BIN_SIZE)]++;
    cumulative[0]=bin[0];
    for(int k=1;k<=200;k++)
        cumulative[k]=cumulative[k-1]+bin[k];
      
    for(int k=0;k<200;k++)
    {
        //printf("%d-%d\t%d\n",k*BIN_SIZE,(k+1)*BIN_SIZE,bin[k]);
        printf("%d\t%.8f\t%.8f\n",k*BIN_SIZE,bin[k]/((cnt+1)*BIN_SIZE),cumulative[k]/(cnt+1));
    }

/*
    for(int k=0;k<=cnt;k++)
    {
        
        if(orientation[k]!=0)
        {
            survival_time_sorted[scnt]=survival_time[k]; //Unsorted survival times
            printf("%d\t%d\t%d\t%d\n",scnt,survival_time[k],survival_time_sorted[k],orientation[k]);
            //printf("%d\n",survival_time[k]);
            scnt++;
        }
        for(int l=0; l < survival_time[k];l++)
        {
            printf("%d\t%d\t%d\t%d\n",time_cnt,orientation[k],survival_time[k],k);
            time_cnt++;
        }
    }
*/
   //}


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

void sort(int n, int arr[])
{
	unsigned long *istack;
        int i,ir=n-1,j,k,l=0;
	int jstack=0;
	int a,temp;

	istack=lvector(1,NSTACK);
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				for (i=j-1;i>=l;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
				}
				arr[i+1]=a;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
                        //printf("Partitioning Element %d\n",k);
                        //printf("SWAP1 arr[%d]=%d arr[%d]=%d\n",k,arr[k],l+1,arr[l+1]);
			SWAP(arr[k],arr[l+1]);
			if (arr[l] > arr[ir]) {
                                //printf("SWAP1 arr[%d]=%d arr[%d]=%d\n",l,arr[l],ir,arr[ir]);
				SWAP(arr[l],arr[ir]);
			}
			if (arr[l+1] > arr[ir]) {
                                //printf("SWAP2 arr[%d]=%d arr[%d]=%d\n",l+1,arr[l+1],ir,arr[ir]);
				SWAP(arr[l+1],arr[ir]);
			}
			if (arr[l] > arr[l+1]) {
                                //printf("SWAP3 arr[%d]=%d arr[%d]=%d\n",l,arr[l],l+1,arr[l+1]);
				SWAP(arr[l],arr[l+1]);
			}
			i=l+1;
			j=ir;
			a=arr[l+1];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
                                //printf("SWAP3 arr[%d]=%d arr[%d]=%d\n",i,arr[i],j,arr[j]);
				SWAP(arr[i],arr[j]);
			}
			arr[l+1]=arr[j];
			arr[j]=a;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_lvector(istack,1,NSTACK);
}

