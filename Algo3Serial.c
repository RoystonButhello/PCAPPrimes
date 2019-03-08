//Algorithm 3.2 - Sequential Prime Sieve Using Fast Marking (2 excluded)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
typedef unsigned __int32 integer;

int main()
{
	//Datatype Range:		[0,4294967296]
	//Acceptable i/p Range: [1,2147483647]
	integer n, bound = (integer)(pow(2, 32) - pow(2, 31) - 1);							
	printf("Enter limit: "); scanf("%d", &n); n--; 
	if (n > bound)
	{
		n = bound;
		printf("Bounded down to: %d\n", n);
	}

	//Boolean array initialized to false
	_Bool *C = (_Bool *)calloc(n, sizeof(_Bool));
	integer i, j;							

	clock_t START_TIME, END_TIME;
	double  CPU_TIME;

	START_TIME = clock();
	for(i=3; i<=(integer)pow(2*n,0.5); i+=2)			//i->[3,5,7,9...,sqrt(2*n)] | i corresponds to C[(i-3)/2]
	{
		if (C[(i-3)/2]==0)		//false implies <i> is prime
		{
			for (j=(integer)pow(i,2); j<=2*n; j+=(2*i))	//j->[i^2,n] | increments by 2*i
				C[(j - 3) / 2] = 1;
			/*LOGIC
			Instead of using modulus operation to check if it is prime within the j-loop,
			change increment from "2" to "2i" and directly mark them as non-prime
			"i^2 + 2k*i" is essentially "m*i", i.e., it is a multiple of i for sure
			*/
		}
	}
	END_TIME = clock();

	/*
	//Print all primes found
	short count = 1;
	printf("Primes in range:\n");
	for(i=0; i<n; i++)
	{
		if (C[i] == 0)	  {printf("%d\t",2*i+3); count++;}
		if (count%10 == 0) { printf("\n"); count /= 10;}
	}
	*/

	/*
	//Print last prime found
	for(i=n-1; i>=0; i--)
		if (C[i] == 0) break;
	printf("Largest prime in range: %d\n",2*i+3);
	*/

	
	short count = 0;
	for (i = n - 1; count < 10; i--)
	{
		if (C[i] == 0) 
		{
			count++;
			printf("%d>%d\n", i, 2 * i + 3);
		}
	}
	

	//Detail execution time
	CPU_TIME = ((double)(END_TIME - START_TIME))/CLOCKS_PER_SEC;
	printf("\nCPU Time: %0.3f seconds\n", CPU_TIME);
	return 0;
}