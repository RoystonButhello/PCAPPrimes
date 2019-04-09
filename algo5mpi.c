//Algorithm 5.1 - Sequential Segmented Prime Sieve Using Fast Marking w/ Block Decomposition
//Phase 1 covers seed values upto sqrt(n)
//Phase 2 covers seed values beyond sqrt(n) in blocks
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>


typedef unsigned long long integer;
const size_t size = sizeof(integer);

//The number of primes to find
integer n;

integer m;
integer blocksize;
integer plen;

//Global declaration of Boolean array
_Bool *mark;
integer *P;

//Log execution time
clock_t START_TIME, END_TIME;
double  CPU_TIME1 = 0.0, CPU_TIME2 = 0.0;

//Utility function to calculate postive integer-powers
integer power(integer val, integer exp)
{
	integer temp = val;
	for (integer i = 1; i < exp; i++) temp *= val;
	return temp;
}

//Utility function to bound seed value to avoid overflow
integer bound(integer n)
{
	short s = size * 8;
	integer limit = power(2, s) - power(2, (s / 2) + 1) - 1;
	if (n > limit)
	{
		n = limit;
		printf("Bounded down to %llu\n", limit);
	}
	return n;
}

//Utility function to approximate no. of primes between 1->n as n/ln(n)
integer trimSize(integer n)
{
	long double e = 2.7183;
	integer exp = 1;
	while (pow(e, exp) < n)
		exp++;
	return n / (exp - 2);
}

//Driver function
int main(int argc, char *argv[])
{
	//Range: Data-type dependent
    	
    int p_rank,p_size;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&p_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&p_size);	
	
    if(p_rank==0)
    {	

       printf("Enter limit: "); 
	   scanf("%llu", &n); 
    
	   n = bound(n);
	   m = sqrt(n);
	   blocksize = m >> 1;
	   plen = trimSize(m);
	
	   if (m % 2 == 0) 
		  m--;

    
       //Boolean array initialized to false
       mark = (_Bool *)calloc(blocksize, sizeof(_Bool));	//Represents [2,3,5,7,9,11,...,sqrt(n)]
	   P = (integer *)calloc(plen, size);
	
	   if (mark == NULL || P == NULL) 
	   { 
		  printf("Memory Allocation Failed!\n"); 
		  exit(1); 
	   }
	

	//Prep file-pointer to write results to text file
	/*FILE *fp = fopen("primes.txt", "w");
	if (fp == NULL) 
	{ 
	printf("File-exception!\n"); 
	exit(0); 
	}*/

	
	
	

	//Loop variables
    integer i,j,k;
	//Setup-Phase: Calculate all primes in the range [3,m]
	START_TIME = clock();
	for (i = 3, k = 0; i < m; i += 2)	//i->[3,5,7,9...,sqrt(n)] | i corresponds to mark[(i-3)/2]
	{
		if (!mark[i >> 1])
		{
			if (i*i <= m)
			{
				for (j = i * i; j <= m; j += (i << 1))	//j->[i^2,n] | increments by 2*i
					mark[j >> 1] = 1;
			}
			printf("%llu\n", i);
			P[k++] = i;
		}
	}
	free(mark);
	plen = k;
	END_TIME = clock();
	CPU_TIME1 = ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;
    
    }
	
	
	//Beginning of segmentation
	//Calculate all primes in the range [m,n] in segmented blocks of size (m/2)
	//Doubled as a blocksize of X covers 2X digits

	//Max limit 1XE10
    
    //Benchmarks 
    //limit tested:1E9
    //Algo5Serial.c:SETUP_PHASE_CPU_TIME:0.021s
    //Algo5Serial.c:COMPUTE_PHASE_CPU_TIME:7.352s
    //algo5mpi.c:SETUP_PHASE_CPU_TIME:0.026s   (after adding if(p_rank==0))
    //algo5mpi.c:COMPUTE_PHASE_CPU_TIME:7.173s (after adding if(p_rank==0))
   
    integer i,j;

	integer min = blocksize << 1;
	integer max = blocksize << 2;
	integer limit = blocksize;
	//Blocks contain ranges [blocksize*k, blocksize*(k+1)]
	START_TIME = clock();
	while (min < n)
	{
		if (max >= n)
		{
			max = n; 
			limit = (max - min) >> 1 + 1;
		}
		
		mark = (_Bool *)calloc(blocksize + 1, sizeof(_Bool));
		
		for (i = 0;i < plen;i++)
		{
			//Find smallest odd-multiple of P[i] w/i range [min,max]
			integer minb = (integer)(min / P[i]) * (P[i]) + P[i];
			if (~minb & 1) minb += P[i];
			//Mark odd-multiples of P[i] as composite
			for (j = minb; j < max; j += 2 * P[i])
				mark[(j - min) >> 1] = 1;
		}
		
		END_TIME = clock();
		CPU_TIME2 += ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;
		
		for (i = 0; i < limit; i++)
			if (!mark[i])
				printf("%llu\n", min + (i << 1) + 1);
		
		START_TIME = clock();
		min += (blocksize << 1);
		max += (blocksize << 1);
		free(mark);
	}
	//fclose(fp);
	
    END_TIME = clock();
	CPU_TIME2 += ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;
	printf("\n\nSETUP-PHASE CPU Time: %0.3f seconds\n", CPU_TIME1);
	printf("COMPUTE-PHASE CPU Time: %0.3f seconds\n", CPU_TIME2);
	
    MPI_Finalize();
	
	return 0;

}