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
	for (;threadCount < 4 ;threadCount++)
	{
		if(total_bytes / threadCount < 8)
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
void sieveBlock(uint32 *, uint64, uint32, uint32, int, _Bool);

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
		bound(&n,&m);
		printf("n = %llu\nm = %llu\n", n, m);
		P = segboolsieve(n, m, &plen);
	}

	//Synchronise processes and broadcast prime-array details
	MPI_Barrier(COMM);
	MPI_Bcast(&plen, 1, MPI_LONG_LONG_INT, 0, COMM);
	MPI_Bcast(&threadCount, 1, MPI_LONG_LONG_INT, 0, COMM);

	//Synchronize processes
	MPI_Barrier(COMM);

	//Space allocated for prime arrays (worst case: 1GB per array allocated)
	if(rank!=0 && rank<threadCount)
		P = (uint32 *)calloc(plen + 1, sizeof(uint32));

	if(rank<threadCount && P == NULL)
	{ printf("\nMemory Allocation Failed!\n"); fflush(stdout); exit(1); }

	//Broadcast primary parameters
	MPI_Bcast(&n, 1, MPI_LONG_LONG_INT, 0, COMM);
	MPI_Bcast(&m, 1, MPI_LONG_LONG_INT, 0, COMM);
	MPI_Bcast(P, plen, MPI_LONG_LONG_INT, 0, COMM);

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
		sieveBlock(P, m, plen, completed, rank, 0);
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
		sieveBlock(P, m, plen, completed, rank, 1);
		END_TIME = clock();
		PROC_TIME += ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;
		fflush(stdout);
		completed = threadCount;
	}

	free(P);
	if(rank<threadCount)
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

void sieveBlock(uint32 *P, uint64 m, uint32 plen, uint32 completed, int rank, _Bool print)
{
	int pid;
	uint32 prime;
	uint64 i, j, minb;
	uint64 min = (uint64)(completed + rank) * m + 1;
	uint64 max = min + m - 2;
	///printf("Segment %3u sieving %llu->%llu | min = (%u + %d)*%llu + 1\n", completed + rank, min, max, completed, rank, m);

	prime = P[0];
	uint64 *mark = (uint64 *)calloc(m / 128, size64);
	for (i = 1;(prime*prime <= max) && (i < plen);i++)
	{
		minb = ((min / prime) * prime);
		if (minb < min) minb += prime;
		if (~minb & 1)	minb += prime;
		for (j = minb;j <= max;j += (prime << 1))
			setBit(mark, j - min + 1);
		prime = P[i];
	}
	if(print)
	for (pid = 0;pid <= rank;pid++)
	{
		if (pid == rank)
		{
			for (j = max;j >= min;j -= 2)
				if (isPrime(mark, j - min + 1))
					break;
			printf("Segment %3u: %llu | %llu\n", completed + rank, j, max);
			break;
		}
	}
	free(mark);
}