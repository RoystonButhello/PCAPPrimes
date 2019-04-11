#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef unsigned __int64 uint64;
typedef unsigned __int32 uint32;

//Utility function to calculate postive uint64-powers
uint64 power(uint64 val, uint64 exp)
{
	uint64 temp = val;
	for (uint64 i = 1; i < exp; i++) temp *= val;
	return temp;
}

//Utility function to approximate no. of primes between 1->n as n/ln(n)
uint32 trimSize(uint64 n)
{
	long double e = 2.7183;
	uint64 exp = 1;
	while (pow(e, exp) < n)
		exp++;
	return (uint32)(n / (exp - 2));
}

//Utility function to round n and m for convenient bit-mapping
void bound(uint64 *n, uint64 *m)
{
	printf("Rounded %llu to ", *n);
	double temp = floor(sqrt(*n));
	temp = ceil(temp / 128) * 128;
	*m = (uint64)temp;
	*n = (*m) * (*m);
	printf("%llu\n", *n);
}

//Utility function to generate an array of primes upto n (worst case: 1GB to store all 32-bit primes)
uint32 *segboolsieve(uint64 n, uint64 m, uint32 *plen)
{
	uint32 blocksize = (uint32)m >> 1;
	*plen = trimSize(m + 1);

	//Boolean array initialized to false
	_Bool *mark = (_Bool *)calloc(blocksize + 1, sizeof(_Bool));	//Represents [2,3,5,7,9,11,...,sqrt(n)]

	//Array to store primes b/w [2,m]
	uint32 *P = (uint32 *)calloc((*plen), sizeof(uint32) + 1);
	if (mark == NULL || P == NULL) { printf("Memory Allocation Failed!\n"); exit(1); }
	uint32 i, j, k, offset;

	//Log execution time
	clock_t START_TIME, END_TIME;
	double  CPU_TIME = 0.0;

	//Setup-Phase: Calculate all primes in the range [3,m]
	START_TIME = clock();
	for (i = 5, k = 1, offset = 2; i < m; i += offset, offset = 6 - offset)	//i->[3,5,7,9...,sqrt(n)] | i corresponds to mark[(i-3)/2]
	{
		if (!mark[i >> 1])
		{
			for (j = i * i; j < m; j += (i << 1))	//j->[i^2,n] | increments by 2*i
				mark[j >> 1] = 1;
			P[k++] = i;
		}
	}
	free(mark);
	P[0] = 3;
	*plen = k;
	END_TIME = clock();
	CPU_TIME = ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;

	printf("\nSETUP-PHASE CPU Time: %0.3f seconds\n", CPU_TIME);
	return P;
}