//mpiexec -n 10 MPITest.exe
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef unsigned __int64 integer;
const size_t size = sizeof(integer);
const integer pow2_32 = 4294967296;
const int threads = 6;

integer power(integer, integer);
integer bound(integer);
_Bool isBitPrime(integer *, integer);
void markBitComposite(integer *, integer);
void sieveBlock(integer *, integer, int);

int main(int argc, char *argv[])
{
	int rank, size, color, temp;
	clock_t START_TIME, END_TIME;
	integer i, j, n = 0, m = 0, plen = 0;
	int smallsieve = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (rank == 0)
	{
		if (size < threads)
		{
			printf("ERROR: Insufficient processes.\n");
			exit(0);
		}
	}
	if (rank < threads)	color = 1;
	else		color = 0;
	MPI_Comm COMM;
	MPI_Comm_split(MPI_COMM_WORLD, color, rank, &COMM);

	if (rank == 0)
	{
		printf("Enter n: ");
		fflush(stdout);
		temp = scanf("%llu", &n);
		n = bound(n);
		m = pow2_32 - 1;	//Use serial sieve for n<2^32
		if (n <= m)
		{
			m = n;
			smallsieve = 1;
			printf("\nUsing small-sieve:-\n\n");
			fflush(stdout);
		}
		else if (n % pow2_32 > 0)	//Round n to nearest multiple of 2^32
		{
			printf("\nRounded %llu to ", n);
			fflush(stdout);
			n = ((n / pow2_32) + 1) * pow2_32;
			printf("%llu\n\n", n);
			fflush(stdout);
		}
		if (~n & 1) n--;
		if (~m & 1) m--;
		plen = (m / 128) + 1;
	}

	//Integer array initialized to 0
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&smallsieve, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&plen, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
	integer *P = (integer *)calloc(plen, size);
	if (P == NULL) { printf("\nMemory Allocation Failed!\n"); fflush(stdout); exit(1); }

	if(rank == 0)
	{
		//Log execution time
		double  CPU_TIME = 0.0;
		printf("BEGINNING SERIAL SIEVE\n");
		fflush(stdout);

		//Setup-Phase: Sieve the range [3,m] w/ a Serial MOD-2 BitSieve
		START_TIME = clock();
		for (i = 3; i*i <= m; i += 2)	//i->[3,5,7,9...,sqrt(n)]
		{
			if (isBitPrime(P, i))
			{
				for (j = i * i; j <= m; j += (i << 1))
					markBitComposite(P, j);
			}
		}
		END_TIME = clock();
		CPU_TIME = ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;

		printf("\nSETUP-PHASE CPU Time: %0.3f seconds\n", CPU_TIME);

		if (smallsieve)
		{
			int count = 0;
			for (i = n; count < 10; i -= 2)
			{
				if (isBitPrime(P, i))
				{
					count++;
					printf("%llu\n", i);
				}
			}
			free(P);
		}
		else printf("SETUP_COMPLETE\n\n");
		fflush(stdout);
	}

	if (smallsieve)
	{
		MPI_Finalize();
		return 0;
	}

	MPI_Bcast(&n, 1, MPI_LONG_LONG_INT, 0, COMM);
	MPI_Bcast(&m, 1, MPI_LONG_LONG_INT, 0, COMM);
	MPI_Bcast(P, (int)plen, MPI_LONG_LONG_INT, 0, COMM);

	MPI_Barrier(MPI_COMM_WORLD);
	double PROC_TIME = 0.0;
	integer max_bytes = pow2_32 / 128;		//Max no. of bytes required to store a segment
	integer len = ((n - m) / 128);			//Total number of bytes to cover ALL values from (pow2_32->n)
	integer segments = len / max_bytes;		//No. of segments
	integer completed = 0;
	integer *BYTES = NULL;

	if (rank == 0)
	{
		printf("SEGMENT PREPERATION COMPLETED.\n");
		printf("Max bytes: \t%llu\n", max_bytes);
		printf("Total bytes: \t%llu\n", len);
		printf("%llu segments\n\n", segments, m, n);
		fflush(stdout);
	}

	MPI_Barrier(COMM);
	while (completed + threads <= segments)
	{
		START_TIME = clock();
		sieveBlock(P, completed, rank);
		END_TIME = clock();
		PROC_TIME += ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;
		MPI_Barrier(COMM);
		completed += threads;
	}
	if (completed < segments)
	{
		integer remaining = segments - completed;
		if (rank < remaining)
		{
			START_TIME = clock();
			sieveBlock(P, completed, rank);
			END_TIME = clock();
			PROC_TIME += ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;
		}
		completed += remaining;
	}
	MPI_Barrier(COMM);
	free(P);
	if(rank<threads)
		printf("\nProcess %d took %f seconds", rank, PROC_TIME);
	MPI_Finalize();
	return 0;
}

//Utility function to calculate postive integer-powers
integer power(integer val, integer exp)
{
	integer temp = val;
	for (integer i = 1; i < exp; i++)
		temp *= val;
	if (temp == 0) temp--;	//Handle overflow to 0
	return temp;
}

//Utility function to bound seed value to avoid overflow
integer bound(integer n)
{
	short s = 64;
	integer limit = power(2, s) - power(2, (s / 2) + 1) - 1;
	if (n > limit)
	{
		n = limit;
		printf("Bounded down to %llu\n", limit);
	}
	return n;
}

_Bool isBitPrime(integer *P, integer x)
{
	return ((P[x / 128] & ((integer)1 << ((x >> 1) & 63)))) ? 0 : 1;
}

//Key-function to set bit corresponding to a number <x>
void markBitComposite(integer *P, integer x)
{
	P[x / 128] |= ((integer)1 << ((x >> 1) & 63));
}

void sieveBlock(integer *P, integer completed, int rank)
{
	integer i, j, minb;
	integer id = completed + rank + 1;
	integer min = (pow2_32 * id) + 1;
	integer max = (pow2_32 * (id + 1)) - 1;
	///printf("Process %llu sieving %llu->%llu\n", id, min, max);
	fflush(stdout);
	integer *BYTES = (integer *)calloc(pow2_32 / 128, size);
	for (i = 3;(i*i) <= max;i += 2)
	{
		if (isBitPrime(P, i)) //If i is prime
		{
			minb = (integer)(min / i)*i + i;
			if (~minb & 1) minb += i;
			//printf("%llu sieving %llu->%llu\n", i,minb,max);
			for (j = minb;j <= max;j += (i << 1))
			{
				markBitComposite(BYTES, j - min + 1);
				//BYTES[temp / 128] |= ((integer)1 << ((temp >> 1) & 63)); //Mark as composite
				//printf("Segment %llu: %llu marked %llu\n", id, i, j);
			}
		}
	}
	FILE *fp = fopen("primes.txt", "w");
	fseek(fp, 0, SEEK_END);
	for (i = 0;i <= rank;i++)
	{
		if (i == rank)
		{
			printf("Last prime in segment %llu: ", id);
			for (j = max;j >= min;j -= 2)
				if (isBitPrime(BYTES, j - min + 1))
					break;
			fprintf(fp,"%llu\n", j);
			printf("%llu\n", j);
			break;
		}
		else
		 fprintf(fp, "\n");
	}
	fclose(fp);
	free(BYTES);
}