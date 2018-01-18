#include "nrutil.h"
#include "nrutil.c"

#define SWAP(a,b); {temp=(a);(a)=(b);(b)=temp;}
#define M 7
#define NSTACK 50

void sort(unsigned long , int []);

int main()
{
    int arr[10]={8,23,89,1,67,45,99,17,1,50};
    sort(10,arr);
    for(int i=0;i<10;i++)
    {
        printf("%d\n",arr[i]);
    }
    printf("%d %d\n",arr[0],arr[1]);
    return 0;
}

void sort(unsigned long n, int arr[])
{
    //for(int i=0;i<n;i++)
      //  printf("%d\t%d\n",i,arr[i]);
	unsigned long *istack;
        int i,ir=n-1,j,k,l=0;
	int jstack=0;
	int a,temp;

	istack=lvector(1,NSTACK);
	for (;;) {
		if (ir-l < M) {
                        printf("Insertion Sort Kicks in l=%d ir=%d\n",l,ir);
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				for (i=j-1;i>=l;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
                                        printf("i=%d l=%d\n",i,l);
				}
				arr[i+1]=a;
                                printf("j=%d i=%d arr[%d]=%d\n",j,i,i+1,arr[i+1]);
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
                        printf("ir=%d l=%d\n",ir,l);
		} else {
			k=(l+ir) >> 1;
                        printf("Partitioning Element %d\n",k);
                        printf("ir=%d l=%d\n",ir,l);
                        printf("SWAP0 arr[%d]=%d arr[%d]=%d\n",k,arr[k],l+1,arr[l+1]);
			SWAP(arr[k],arr[l+1]);
			if (arr[l] > arr[ir]) {
                                printf("SWAP1 arr[%d]=%d arr[%d]=%d\n",l,arr[l],ir,arr[ir]);
				SWAP(arr[l],arr[ir]);
			}
			if (arr[l+1] > arr[ir]) {
                                printf("SWAP2 arr[%d]=%d arr[%d]=%d\n",l+1,arr[l+1],ir,arr[ir]);
				SWAP(arr[l+1],arr[ir]);
			}
			if (arr[l] > arr[l+1]) {
                                printf("SWAP3 arr[%d]=%d arr[%d]=%d\n",l,arr[l],l+1,arr[l+1]);
				SWAP(arr[l],arr[l+1]);
			}
			i=l+1;
			j=ir;
			a=arr[l+1];//Partitioning Element
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
                                printf("SWAPCROSS arr[%d]=%d arr[%d]=%d\n",i,arr[i],j,arr[j]);
				SWAP(arr[i],arr[j]);
			}
                        printf("PE=%d inserted now arr[%d]=%d arr[%d]=%d \n",a,l+1,arr[l+1],j,arr[j]);
			arr[l+1]=arr[j];//Insert partitioning element
			arr[j]=a;
			jstack += 2;//Number of sub arrays 
                        printf("jstack=%d ir=%d i=%d j=%d l=%d\n",jstack,ir,i,j,l);
			if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
                                printf("Case1 ir-i+1>= j-l: istack[%d]=%d istack[%d]=%d ir=%d istack[jstack]=%d istack[jstack-1]=%d",jstack,ir,jstack-1,i,ir,istack[jstack],istack[jstack-1]);
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
                                printf("Case2 ir-i+1 < j-l: istack[%d]=%d istack[%d]=%d l=%d istack[jstack]=%d istack[jstack-1]=%d\n",jstack,j-1,jstack-1,l,i,istack[jstack],istack[jstack-1]);
			}
		}
	}
	free_lvector(istack,1,NSTACK);
}


