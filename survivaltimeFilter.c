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
#define MAXLENGTH 1000000
#define PERIOD 10000
#define RUNMAX 201
#define NMAX 2000
#define MAXFRAMES 20001
#define M 7
#define NSTACK 60

double cnode[RUNMAX][MAXFRAMES];
int survival_time[MAXLENGTH],orientation[MAXLENGTH];
int filter_survival_time[MAXLENGTH],filter_orientation[MAXLENGTH];
int survival_time_sorted[MAXLENGTH],survival_time_unsorted[MAXLENGTH];

int NX,NY,RUNS,STEPS,FRAMES;
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
   case 7:
       NX = atoi(argv[1]);
       NY = atoi(argv[2]);
       KAPPA = atof(argv[3]);
       STEPS = atoi(argv[4]);
       FILTER1 = atof(argv[5]);
       FILTER2 = atof(argv[6]);
       break;
   default:
       print_and_exit("Usage Pass command line arguments:NX NY Kappa STEPS FILTER1 FILTER2\n");
   }

   FRAMES = STEPS/PERIOD;

   sprintf(data_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/cnode.bin",NX,NY,KAPPA);
   printf("%s\n",data_file);
   if(NULL==(Fin=fopen(data_file,"rb")))
   	print_and_exit("I could not open binary file %s\n",data_file);
   

   //printf("frames %d\n",frames);

   //Initializing the output file
   sprintf(out_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/survivaltime1.log",NX,NY,KAPPA);
   printf("Output file: %s\n",out_file);
   if(NULL==(Fout=fopen(out_file,"w")))
         print_and_exit("I could not open file with simulation survival time %s\n",out_file);

   //We read the data file
   fread(&RUNNUM,sizeof(int),1,Fin);
   printf("Total Runs: %d\n",RUNNUM);
   RUNNUM=1;//reading just the 1st run for checking 
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

        //Converting to Up, Down and Transition
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
    

    // Eliminating small excursions to Transition from Up/Down uses FILTER1
    for (int k=1;k<cnt;k++)
    {
      if(((orientation[k-1]==-1 && orientation[k+1]==-1 && survival_time[k-1]>10) || (orientation[k-1]==1 && orientation[k+1]==1 && survival_time[k-1]>10)) && orientation[k]==0)
      {
          //printf("Excursion Survival time %d : %d\n",i,survival_time[k]);
          if(survival_time[k]<FILTER1)
          {
              survival_time[k-1]+=survival_time[k]+survival_time[k+1];
              left_shift_array(survival_time,cnt-2,k);
              left_shift_array(orientation,cnt-2,k);
              cnt-=2;//since we have shifted two places left
              k--;//Next k value remains same due to shift
          }
      }
    }

    // Eliminating small excursions to Up/Down from Transition used FILTER2
    for (int k=1;k<cnt;k++)
    {
        if(orientation[k-1]==0 && orientation[k+1]==0 && (orientation[k]==1  || orientation[k]==-1) && survival_time[k-1]>10)
        {
            //printf("Survival time %d : %d at %d orientation %d\n",k,survival_time[k],j,orientation[k]);
            if(survival_time[k]<FILTER2)
            {
              survival_time[k-1]+=survival_time[k]+survival_time[k+1];
              left_shift_array(survival_time,cnt-2,k);
              left_shift_array(orientation,cnt-2,k);
              cnt-=2;//since we have shifted two places left
              k--;//Next k value remains same due to shift
            }
        }
    }


    // Eliminating small excursions to Transition from Up/Down uses FILTER1
    for (int k=1;k<cnt;k++)
    {
      if(((orientation[k-1]==-1 && orientation[k+1]==-1 && survival_time[k-1]>10) || (orientation[k-1]==1 && orientation[k+1]==1 && survival_time[k-1]>10)) && orientation[k]==0)
      {
          //printf("Excursion Survival time %d : %d\n",i,survival_time[k]);
          if(survival_time[k]<FILTER1)
          {
              survival_time[k-1]+=survival_time[k]+survival_time[k+1];
              left_shift_array(survival_time,cnt-2,k);
              left_shift_array(orientation,cnt-2,k);
              cnt-=2;//since we have shifted two places left
              k--;//Next k value remains same due to shift
          }
      }
    }

    // Eliminating small excursions to Up/Down from Transition uses FILTER2
    for (int k=1;k<cnt;k++)
    {
        if(orientation[k-1]==0 && orientation[k+1]==0 && (orientation[k]==1  || orientation[k]==-1) && survival_time[k-1]>10)
        {
            //printf("Survival time %d : %d at %d orientation %d\n",k,survival_time[k],j,orientation[k]);
            if(survival_time[k]<FILTER2)
            {
              survival_time[k-1]+=survival_time[k]+survival_time[k+1];
              left_shift_array(survival_time,cnt-2,k);
              left_shift_array(orientation,cnt-2,k);
              cnt-=2;//since we have shifted two places left
              k--;//Next k value remains same due to shift
            }
        }
    }
   }

    //Printing orientation after passing through the filters 
    int time_cnt=0.2*FRAMES;
    for(int k=0;k<=cnt;k++)
    {
       for(int l=0; l < survival_time[k];l++)
        {
            printf("%d\t%d\t%d\t%d\n",time_cnt,orientation[k],survival_time[k],k);
            time_cnt++;
        }
    }
 




/*
    //printf("cnt %d\n",cnt);
    int time_cnt=0.2*FRAMES;
    for(int k=0;k<=cnt;k++)
    {
        if(orientation[k]!=0)
        {
            survival_time_unsorted[scnt]=survival_time[k]; //Unsorted survival times
            survival_time_sorted[scnt]=survival_time_unsorted[scnt];//This array will be sorted
            //printf("%d\t%d\t%d\n",scnt,survival_time_unsorted[scnt],survival_time_sorted[scnt]);
            scnt++;
        }
    }
   }


    //Sort the survival times
    sort(scnt,survival_time_sorted);
    printf("scnt = %d\n",scnt);    
    for(int k=0;k<scnt;k++)
    {
        printf("%d\t%d\t%d\n",k,survival_time_unsorted[k],survival_time_sorted[k]);
    }

    //Bin the sorted survival times
    int BIN_SIZE=5;
    int bin[MAXLENGTH];
    int bin_cnt=1,bin_entry=0;
    for(int k=0;k<scnt;k++)
    {
        if (survival_time_sorted[k] <= bin_cnt*BIN_SIZE)
        {
            bin_entry++;
            //printf("bin_cnt = %d\n",bin_cnt);
            if(k==scnt-1)
                bin[bin_cnt-1]=bin_entry;
        }
        else
        {
            bin[bin_cnt-1]=bin_entry;
            bin_cnt++;
            bin_entry=0;
            k--;
        }
    }

    for(int k=0;k<bin_cnt;k++)
    {
        //printf("%d-%d\t%d\n",k*BIN_SIZE,(k+1)*BIN_SIZE,bin[k]);
        //printf("%d\t%d\n",k,bin[k]);
    }
   

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

