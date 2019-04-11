#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
typedef unsigned __int64 integer;
const size_t size = sizeof(integer);

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

integer *segboolsieve(integer n, integer m, integer &plen, bool smallsieve)
{
	integer blocksize = m >> 1;
	if (smallsieve)
		plen = trimSize(n);
	else
		plen = trimSize(m + 2);

	//Boolean array initialized to false
	bool *mark = (bool *)calloc(blocksize + 1, sizeof(bool));	//Represents [2,3,5,7,9,11,...,sqrt(n)]

	//Array to store primes b/w [2,m]
	integer *P = (integer *)calloc(plen, (size_t)size + 1);
	if (mark == NULL || P == NULL) { printf("Memory Allocation Failed!\n"); exit(1); }
	integer i, j, k, offset;

	//Log execution time
	clock_t START_TIME, END_TIME;
	double  CPU_TIME1 = 0.0, CPU_TIME2 = 0.0;

	//Setup-Phase: Calculate all primes in the range [3,m]
	START_TIME = clock();
	for (i = 5, k = 1, offset = 2; i <= m; i += offset, offset = 6 - offset)	//i->[3,5,7,9...,sqrt(n)] | i corresponds to mark[(i-3)/2]
	{
		if (!mark[i >> 1])
		{
			if (i*i <= n)
			{
				for (j = i * i; j <= m; j += (i << 1))	//j->[i^2,n] | increments by 2*i
					mark[j >> 1] = true;
			}
			P[k++] = i;
		}
	}
	free(mark);
	P[0] = 3;
	plen = k;
	END_TIME = clock();
	CPU_TIME1 = ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;

	//If CUDA Sieve will be used
	if (!smallsieve)
	{
		printf("\nSETUP-PHASE CPU Time: %0.3f seconds\n", CPU_TIME1);
		return P;
	}

	//Calculate all primes in the range [m,n] in segmented blocks of size (m/2)
	//Doubled as a blocksize of X covers 2X digits
	integer min = (blocksize << 1) + 1;
	integer max = (blocksize << 2) - 1;
	//Blocks contain ranges [blocksize*k, blocksize*(k+1)]
	START_TIME = clock();
	while (min < n)
	{
		if (max > n) max = n;
		mark = (bool *)calloc(blocksize + 1, sizeof(bool));
		for (i = 0;i < plen;i++)
		{
			//Find smallest odd-multiple of P[i] w/i range [min,max]
			integer minb = (integer)(min / P[i]) * (P[i]) + P[i];
			if (~minb & 1) minb += P[i];
			//Mark odd-multiples of P[i] as composite
			for (j = minb; j <= max; j += 2 * P[i])
				mark[(j - min) >> 1] = true;
		}
		for (i = min; i <= max;i += 2)
		{
			if (!mark[(i - min) >> 1])
				P[k++] = i;
		}
		min += (blocksize << 1);
		max += (blocksize << 1);
		free(mark);
	}
	plen = k;
	END_TIME = clock();
	CPU_TIME2 += ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;

	printf("\nSETUP-PHASE CPU Time: %0.3f seconds\n", CPU_TIME1 + CPU_TIME2);
	return P;
}