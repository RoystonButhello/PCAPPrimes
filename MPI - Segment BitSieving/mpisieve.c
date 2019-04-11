//mpiexec -n 4 MPITest.exe
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "smallsieve.h"

typedef unsigned __int64 uint64;
typedef unsigned __int32 uint32;
const size_t size64 = sizeof(uint64);

//Utility function to calculate neccessary number of threads
uint32 threadCalc(uint64 m)
{
	uint64 total_bytes = (m / 128) * 8;
	int threadCount = 1;
	for (;threadCount < 4;threadCount++)
	{
		if (total_bytes / threadCount < 8)
		{
			threadCount--;
			break;
		}
	}
	printf("Using %u threads [%llu bytes per thread]\n", threadCount, total_bytes / threadCount);
	return threadCount;
}

_Bool isPrime(uint64 *, uint64);
void setBit(uint64 *, uint64);
void sieveBlock(uint32 *, uint64, uint32, uint32, int, int);

int main(int argc, char *argv[])
{
	int rank, size, color, temp, threadCount;
	clock_t START_TIME, END_TIME;
	uint64 n = 0, m = 0;
	uint32 plen = 0, *P = NULL;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	threadCount = size;

	if (rank < threadCount)	color = 1;
	else color = 0;
	MPI_Comm COMM;
	MPI_Comm_split(MPI_COMM_WORLD, color, rank, &COMM);

	//Read input from user
	if (rank == 0)
	{
		printf("Enter n: ");
		fflush(stdout);
		temp = scanf("%llu", &n);
		bound(&n, &m);
		printf("n = %llu\nm = %llu\n", n, m);
		P = segboolsieve(n, m, &plen);
	}

	//Synchronise processes and broadcast prime-array details
	MPI_Barrier(COMM);
	MPI_Bcast(&plen, 1, MPI_UNSIGNED_LONG, 0, COMM);
	MPI_Bcast(&threadCount, 1, MPI_INT, 0, COMM);

	//Synchronize processes
	MPI_Barrier(COMM);

	//Space allocated for prime arrays (worst case: 1GB per array allocated)
	if (rank != 0 && rank < threadCount)
		P = (uint32 *)calloc(plen + 1, sizeof(uint32));

	if (rank < threadCount && P == NULL)
	{
		printf("\nMemory Allocation Failed!\n"); fflush(stdout); exit(1);
	}

	//Broadcast primary parameters
	MPI_Bcast(&n, 1, MPI_LONG_LONG_INT, 0, COMM);
	MPI_Bcast(&m, 1, MPI_LONG_LONG_INT, 0, COMM);
	MPI_Bcast(P, plen, MPI_UNSIGNED_LONG, 0, COMM);

	//Initialize secondary parameters
	double PROC_TIME = 0.0;
	uint32 completed = 1;
	uint32 segments = (uint32)(n / m); //segments = 1 for all extraneous processes

	//Print segment details
	if (rank == 0)
	{
		printf("\n%u segments(s) for [%llu->%llu]\n", segments - 1, m + 1, n);
		fflush(stdout);
	}

	//Synchronize processes
	MPI_Barrier(COMM);

	//Segment sieving begins
	while (completed + threadCount <= segments)
	{
		START_TIME = clock();
		sieveBlock(P, m, plen, completed, rank, threadCount);
		END_TIME = clock();
		PROC_TIME += ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;
		fflush(stdout);
		completed += threadCount;
		/*if (rank == 0)
		{
			printf("Segment %u done\n", completed);
			fflush(stdout);
		}*/
	}
	if (completed + rank + 1 <= segments)
	{
		START_TIME = clock();
		sieveBlock(P, m, plen, completed, rank, threadCount);
		END_TIME = clock();
		PROC_TIME += ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;
		fflush(stdout);
		completed = threadCount;
	}

	free(P);
	if (rank < threadCount)
		printf("\nProcess %d took %f seconds", rank, PROC_TIME);

	MPI_Finalize();
	return 0;
}

_Bool isPrime(uint64 *P, uint64 x)
{
	return ((P[x / 128] & ((uint64)1 << ((x >> 1) & 63)))) ? 0 : 1;
}

//Key-function to set bit corresponding to a number <x>
void setBit(uint64 *P, uint64 x)
{
	P[x / 128] |= ((uint64)1 << ((x >> 1) & 63));
}

void sieveBlock(uint32 *P, uint64 m, uint32 plen, uint32 completed, int rank, int threadCount)
{
	int pid;
	uint32 prime;
	uint64 i, j, minb;
	uint64 min = (uint64)(completed + rank) * m + 1;
	uint64 max = min + m - 2;

	prime = P[0];
	uint64 *mark = (uint64 *)calloc(m >> 7, size64);
	for (i = 1;(i < plen) & (prime*prime <= max);i++)
	{
		minb = ((min / prime) * prime);
		if (minb < min) minb += prime;
		if (~minb & 1)	minb += prime;
		for (j = minb;j <= max;j += (prime << 1))
			setBit(mark, j - min + 1);
		prime = P[i];
	}
	if (max > (m*(m - threadCount)))
	{
		for (j = max;j >= min;j -= 2)
			if (isPrime(mark, j - min + 1))
				break;
		printf("%d: Segment %3u: %llu | %llu\n", rank, completed + rank, j, max);
	}
	free(mark);
}

//BENCHMARKS
// SERIAL SEGMENTED BOOLSIEVE
//	n		RAM	    Time
// E08	   620KB   0.206s
// E09	   704KB   1.895s
// E10	   668KB   20.02s 
// E11     904KB   205.2s

//MPI SEGMENTED BITSIEVE (Time taken across all processes averaged
//	n	Processes   AvgTime	
// E08	    1       0.188s
// E08	    2       0.096s
// E08	    4		0.052s
// E08		8		0.036s
// E09	    1       1.740s
// E09	    2       0.892s
// E09	    4		0.532s
// E09		8		0.334s
// E10	    1       17.55s
// E10	    2        9.03s
// E10	    4		 5.18s
// E10		8		 3.38s
// E10		16		3.342s //Point of diminishing returns for 4/8 i7 processor
// E11	    4		 5.18s
// E11		8		 3.38s

