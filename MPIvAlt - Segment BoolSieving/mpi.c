//mpiexec -n 4 MPITest.exe
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

typedef unsigned __int64 uint64;
typedef unsigned __int32 uint32;
const size_t size64 = sizeof(uint64);

//Benchmarks 
//limit tested:1E9
//Algo5Serial.c:SETUP_PHASE_CPU_TIME:0.021s
//Algo5Serial.c:COMPUTE_PHASE_CPU_TIME:7.352s
//algo5mpi.c:SETUP_PHASE_CPU_TIME:0.026s   (after adding if(rank==0))
//algo5mpi.c:COMPUTE_PHASE_CPU_TIME:7.173s (after adding if(rank==0))

//limit tested:1E10
//Algo5Serial.c:SETUP_PHASE_CPU_TIME:0.025s
//Algo5Serial.c:COMPUTE_PHASE_CPU_TIME:65.218s
//algo5mpi.c:SETUP_PHASE_CPU_TIME:0.025s
//algo5mpi.c:COMPUTE_PHASE_CPU_TIME:67.384s

//limit tested:1E6
//algo5mpi.c:SETUP_PHASE_CPU_TIME:0.002s without MPI_Bcast and with Bcast 
//algo5mpi.c:COMPUTE_PHASE_CPU_TIME:0.009s without MPI_Bcast and with Bcast

//limit tested:1E8
//algo5mpi.c:SETUP_PHASE_CPU_TIME:0.006s without MPI_Bcast  
//algo5mpi.c:COMPUTE_PHASE_CPU_TIME:0.772s without MPI_Bcast 
//algo5mpi.c:SETUP_PHASE_CPU_TIME:0.005s with MPI_Bcast  
//algo5mpi.c:COMPUTE_PHASE_CPU_TIME:0.717s with MPI_Bcast 


//The number of primes to find
uint64 n, m;
uint32 blocksize, plen;

//Log execution time
clock_t START_TIME, END_TIME;
double  CPU_TIME1 = 0.0, CPU_TIME2 = 0.0;

//Utility function to calculate postive uint64-powers
uint64 power(uint64 val, uint64 exp)
{
	uint64 temp = val;
	for (uint64 i = 1; i < exp; i++) temp *= val;
	return temp;
}

//Utility function to bound seed value to avoid overflow
void bound(uint64 *n, uint64 *m)
{
	printf("Rounded %llu to ", *n);
	*m = (uint64)(floor(sqrt(*n)));
	*n = (*m) * (*m);
	printf("%llu\n", *n);
	printf("sqrt(%llu) = %llu\n", *n, *m);
}

//Utility function to approximate no. of primes between 1->n as n/ln(n)
uint32 trimSize(uint64 n)
{
	long double e = 2.7183;
	uint64 exp = 1;
	while (pow(e, exp) < n)
		exp++;
	return n / (exp - 2);
}

//Driver function
int main(int argc, char *argv[])
{
	//Range: Data-type dependent

	uint32 *P = NULL;
	_Bool *mark;
	int rank, threadCount;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &threadCount);

	if (rank == 0)
	{
		printf("Enter limit: ");
		fflush(stdout);
		scanf("%llu", &n);

		bound(&n, &m);
		blocksize = m >> 1;
		plen = trimSize(m);
		
		//Boolean array initialized to false
		mark = (_Bool *)calloc(blocksize, sizeof(_Bool));	//Represents [2,3,5,7,9,11,...,sqrt(n)]
		P = (uint32 *)calloc(plen, sizeof(uint32));

		if (mark == NULL || P == NULL)
		{
			printf("Memory Allocation Failed!\n");
			exit(1);
		}

		//Loop variables
		uint64 i, j, k;
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
				P[k++] = i;
			}
		}
		free(mark);
		plen = k;
		END_TIME = clock();
		CPU_TIME1 = ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;
		fflush(stdout);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&plen, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&blocksize, 1, MPI_LONG_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank != 0 && rank < threadCount)
		P = (uint32 *)calloc(trimSize(m), sizeof(uint32));

	if (rank < threadCount && P == NULL)
	{
		printf("\nMemory Allocation Failed!\n"); fflush(stdout); exit(1);
	}
	else printf("P%d allocated memory\n", rank);

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(P, plen, MPI_LONG_INT, 0, MPI_COMM_WORLD);
	for (int i = 0;i < 5;i++)
		printf("P%d: %u\n", rank, P[i]);
	fflush(stdout);

	MPI_Barrier(MPI_COMM_WORLD);
	//Calculate all primes in the range [m,n] in segmented blocks
	uint64 i, j;
	uint32 completed = 1;
	uint32 segments = (uint32)(m);
	uint64 min = (uint64)(completed + rank) * m + 1;
	uint64 max = min + m - 2;
	uint64 limit = blocksize;
	uint64 boundary = segments - 2;
	_Bool print = 0;
	if (~n & 1) n--;

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0)
	{
		printf("\n%u segments(s) for [%llu->%llu]\n", segments - 1, m + 1 , n);
		fflush(stdout);
	}

	START_TIME = clock();
	while (completed + threadCount <= segments)
	{
		//printf("Segment %d at [%10llu,%10llu]\n", completed + rank, min, max);
		fflush(stdout);
		if (completed + threadCount > boundary)
		{
			max = n;
			limit = ((max - min) >> 1) + 1;
			print = 1;
			fflush(stdout);
		}

		mark = (_Bool *)calloc(blocksize + 1, sizeof(_Bool));
		
		for (i = 0;i < plen;i++)
		{
			uint64 minb = (uint64)(min / P[i]) * (P[i]) + P[i];
			if (~minb & 1) minb += P[i];
			for (j = minb; j < max; j += 2 * P[i])
				mark[(j - min) >> 1] = 1;
		}
		END_TIME = clock();
		CPU_TIME2 += ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;

		if (print)
			for (i = 0; i < limit; i++)
				if (!mark[i])
					printf("P%d: %llu\n", rank, min + (i << 1));

		START_TIME = clock();
		min = (uint64)(completed + rank) * m + 1;
		max = min + m - 2;
		completed += threadCount;
		free(mark);
	}
	END_TIME = clock();
	CPU_TIME2 += ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;
	printf("\n\nP%d SETUP-PHASE CPU Time: %0.3f seconds\n", rank, CPU_TIME1);
	printf("P%d COMPUTE-PHASE CPU Time: %0.3f seconds\n", rank, CPU_TIME2);

	MPI_Finalize();
	return 0;

}