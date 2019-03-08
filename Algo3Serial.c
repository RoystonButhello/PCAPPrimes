//Algorithm 3.2 - Sequential Prime Sieve Using Fast Marking (2 excluded)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
typedef unsigned __int32 integer;

integer power(integer val, integer exp)
{
	integer temp = val;
	for (integer i = 1; i < exp; i++) temp *= val;
	return temp;
}

integer bound(integer n)
{
	short size = sizeof(integer) * 8 - 1;
	integer limit = power(2, size) - power(2, (size / 2) + 1) - 1;
	if (n > limit)
	{
		n = limit;
		printf("Bounded down to %u\n", limit);
	}
	return n;
}

int main()
{
	//Datatype Range:		[0,4294967296]
	//Acceptable i/p Range: [1,2147483648]	
	integer n; printf("Enter limit: "); scanf("%u", &n);
	n = bound(n);


	//Boolean array initialized to false
	_Bool *C = (_Bool *)calloc(n, sizeof(_Bool));
	if (C == NULL) { printf("Memory Allocation Failed!\n"); exit(1); }
	integer i, j;

	//Log execution time
	clock_t START_TIME, END_TIME;
	double  CPU_TIME;

	START_TIME = clock();
	for (i = 3; i <= (integer)pow(2 * n, 0.5); i += 2)					//i->[3,5,7,9...,sqrt(2*n)] | i corresponds to C[(i-3)/2]
	{
		if (C[(i - 3) / 2] == 0)		//false implies <i> is prime
		{
			for (j = (integer)pow(i, 2); j <= 2 * n; j += (2 * i))	//j->[i^2,2*n] | increments by 2*i
				C[(j - 3) / 2] = 1;
			/*LOGIC
			Instead of using modulus operation to check if it is prime within the j-loop,
			change increment from "2" to "2i" and directly mark them as non-prime
			"i^2 + 2k*i" is essentially "m*i", i.e., it is a multiple of i for sure
			*/
			//printf("%u stopped at %u\n", i, j);
		}
	}
	END_TIME = clock();

	/*
	//Print all primes found
	short count = 1;
	printf("\nPrimes in range:\n");
	for(i=0; i<n; i++)
	{
		if (C[i] == 0)	  {printf("%u\t",2*i+3); count++;}
		if (count%10 == 0) { printf("\n"); count /= 10;}
	}
	*/

	/*
	//Print last prime found
	for(i=n-1; i>=0; i--)
		if (C[i] == 0) break;
	printf("Largest prime in range: %u\n",2*i+3);
	*/


	printf("\nLast 10 primes:\n");
	short count = 0;
	for (i = n - 1; count < 10; i--)
	{
		if (C[i] == 0)
		{
			count++;
			printf("%u\n", 2 * i + 3);
		}
	}

	//Detail execution time
	CPU_TIME = ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;
	printf("\nCPU Time: %0.3f seconds\n", CPU_TIME);
	return 0;
}