//Algorithm 3.2 - Sequential Prime Sieve Using Fast Marking (2 excluded)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
typedef unsigned __int64 integer;

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
	short size = sizeof(integer) * 8;
	integer limit = power(2, size) - power(2, (size / 2) + 1) - 1;
	if (n > limit)
	{
		n = limit;
		printf("Bounded down to %llu\n", limit);
	}
	return n;
}

//Statistics:
//	   BOOLEAN SIEVE
//	n		RAM		Time
// E07	   520KB   .012s
// E08		48MB   .431s
// E09	   478MB   5.47s
// E10	   4.7GB   74.3s
// E11     calloc failed

//Stats logged via Visual Studio Performance Profiler on i7 4790K @4.00GHz w/ 16GB DDR3 RAM
//Can't take n>E10 due to time and space inefficiency

//Driver function
int main()
{
	//Range: [0,4294836223]
	integer n; printf("Enter limit: "); scanf("%llu", &n);
	n = bound(n);
	integer len = n>>1;

	//Boolean array initialized to false
	_Bool *C = (_Bool *)calloc(len, sizeof(_Bool));	//Represents [2,3,5,7,9,11,...]
	if (C == NULL) { printf("Memory Allocation Failed!\n"); exit(1); }
	integer i, j;

	//Log execution time
	clock_t START_TIME, END_TIME;
	double  CPU_TIME;

	START_TIME = clock();
	for (i = 3; i*i <= n; i += 2)	//i->[3,5,7,9...,sqrt(n)] | i corresponds to C[(i-3)/2]
	{
		if (C[i/2] == 0)			//false implies <i> is prime
		{
			for (j = i * i; j <= n; j += (i << 1))	//j->[i^2,n] | increments by 2*i
			{
				C[j / 2] = 1;
				//printf("%llu marked off\n", j);
			}
			/*LOGIC
			Instead of using modulus operation to check if it is prime within the j-loop,
			change increment from "2" to "2i" and directly mark them as non-prime
			"i^2 + 2k*i" is essentially "m*i", i.e., it is a multiple of i for sure
			*/
			if(i>99000) printf("%llu stopped at %llu\n", i, j - (i << 1));
		}
	}
	END_TIME = clock();
	/*
	//Print all primes found
	short count = 1;
	printf("\nPrimes in range:\n2\t");
	for(i=1; i<len; i++)
	{
		if (C[i] == 0)	  {printf("%llu\t",2*i+1); count++;}
		if (count%10 == 0) { printf("\n"); count /= 10;}
	}
	*/
	
	printf("\nLast 10 primes:\n");
	short count = 0;
	for (i = len - 1; count < 10; i--)
	{
		if (C[i] == 0)
		{
			count++;
			printf("%llu\n", 2 * i + 1);
		}
	}
	
	/*
	integer pcount = 0;
	for (i = 0, pcount = 0; i < len; i++)
		if (C[i] == 0) pcount++;
	printf("%llu primes found\n", pcount);
	*/
	
	//Detail execution time
	CPU_TIME = ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;
	printf("\nCPU Time: %0.3f seconds\n", CPU_TIME);
	free(C);
	return 0;
}
