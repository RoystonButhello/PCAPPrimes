	//Algorithm 5.2 - Sequential Segmented Prime Sieve w/ Block Decomposition and Mod-6 Wheel for setup phase
//Phase 1 covers seed values upto sqrt(n)
//Phase 2 covers seed values beyond sqrt(n) in blocks
//Writing results to file takes significantly longer than both phases combined; selectively write results as needed
//Performance improves over repeated tests due to cache-concurrency
#include <stdio.h>
#include <string.h>
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

//Statistics:
//	     BOOLSIEVE
//	n		RAM		Time
// E07	   520KB   .012s
// E08		48MB   .431s
// E09	   478MB   5.47s
// E10	   4.7GB   74.3s
// E11     calloc failed

//	     BITSIEVE
//	n		RAM		Time
// E07	   532KB   .009s
// E08	   6.5MB   .129s
// E09		60MB   3.35s
// E10	   597MB  44.02s
// E11	   5.8GB     ?
// E12     calloc failed

//  SEGMENTED BOOLSIEVE
//	n		RAM	    Time
// E07	   552KB   0.026s
// E08	   620KB   0.206s
// E09	   704KB   1.895s
// E10	   668KB   20.02s 
// E11     904KB   205.2s

//Stats logged via Visual Studio Performance Profiler on i7 4790K @4.00GHz w/ 16GB DDR3 RAM
//Can't take n>E11 due to time inefficiency; to be further improved via paralellization

//Driver function
int main()
{
	//Range: Data-type dependent
	integer n; printf("Enter limit: "); scanf("%llu", &n); n = bound(n);
	integer m = sqrt(n);
	integer blocksize = m >> 1;
	integer plen = trimSize(m);
	if (m % 2 == 0) m--;

	//Prep file-pointer to write results to text file
	FILE *fp = fopen("primes.txt", "w");
	if (fp == NULL) { printf("File-exception!\n"); exit(0); }

	//Boolean array initialized to false
	_Bool *mark = (_Bool *)calloc(blocksize + 1, sizeof(_Bool));	//Represents [2,3,5,7,9,11,...,sqrt(n)]

	//Array to store primes b/w [2,m]
	integer *P = (integer *)calloc(plen, (size_t)size + 1);
	if (mark == NULL || P == NULL) { printf("Memory Allocation Failed!\n"); exit(1); }
	integer i, j, k, offset;

	//Log execution time
	clock_t START_TIME, END_TIME, F_START_TIME, F_END_TIME;
	double  CPU_TIME1 = 0.0, CPU_TIME2 = 0.0, F_CPU_TIME = 0.0;;

	//Setup-Phase: Calculate all primes in the range [3,m]
	START_TIME = clock();
	for (i = 5, k = 1, offset = 2; i < m; i += offset, offset = 6 - offset)	//i->[3,5,7,9...,sqrt(n)] | i corresponds to mark[(i-3)/2]
	{
		if (!mark[i >> 1])
		{
			if (i*i <= n)
			{
				for (j = i * i; j <= m; j += (i << 1))	//j->[i^2,n] | increments by 2*i
					mark[j >> 1] = 1;
			}
			/*
			F_START_TIME = clock();
			fprintf(fp, "%llu\n", i);
			F_END_TIME = clock();
			F_CPU_TIME += ((double)(F_END_TIME - F_START_TIME)) / CLOCKS_PER_SEC;
			*/
			P[k++] = i;
		}
	}
	free(mark);
	P[0] = 3;
	plen = k;
	END_TIME = clock();
	CPU_TIME1 = ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;

	//Calculate all primes in the range [m,n] in segmented blocks of size (m/2)
	//Doubled as a blocksize of X covers 2X digits
	integer min = blocksize << 1;
	integer max = blocksize << 2;
	integer limit = blocksize;
	//Blocks contain ranges [blocksize*k, blocksize*(k+1)]
	START_TIME = clock();
	while (min < n)
	{
		if (max >= n)
		{
			max = n; limit = (max - min) >> 1 + 1;
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
		//END_TIME = clock();
		//CPU_TIME2 += ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;
		/*
		F_START_TIME = clock();*/
		if(max==n)
			for (i = 0; i < limit; i++)
				if (!mark[i])
					fprintf(fp, "%llu\n", min + (i << 1) + 1);
		/*F_END_TIME = clock();
		F_CPU_TIME += ((double)(F_END_TIME - F_START_TIME)) / CLOCKS_PER_SEC;
		*/
		//START_TIME = clock();
		min += (blocksize << 1);
		max += (blocksize << 1);
		free(mark);
	}
	fclose(fp);
	free(P);
	END_TIME = clock();
	CPU_TIME2 += ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;
	printf("\n\nSETUP-PHASE CPU Time: %0.3f seconds\n", CPU_TIME1);
	printf("COMPUTE-PHASE CPU Time: %0.3f seconds\n", CPU_TIME2);
	printf("FILE_WRITE CPU Time: %0.3f seconds\n", F_CPU_TIME);
	return 0;

}