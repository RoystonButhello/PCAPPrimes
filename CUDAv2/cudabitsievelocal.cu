#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "smallsieve.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef unsigned __int64 uint64;
typedef unsigned __int32 uint32;

///////////////////////////KERNELS START///////////////////////////
//Checks if x is prime; if bit corresponding to x is 0, then return true.
__device__ bool isPrime(uint64 *mark, uint64 x)
{
	return (mark[x / 128] & ((uint64)1 << ((x >> 1) & 63))) ? 0 : 1;
}

//Set the bit corresponding to x
__device__ void bitSet(uint64 *mark, uint64 x)
{
	mark[x / 128] |= ((uint64)1 << ((x >> 1) & 63));
}

__global__ void SieveBlock(uint32 *P, uint64 *mark, uint64 m, uint32 completed, uint32 plen, uint32 seglen)
{
	//Each thread sieves [(pow2_32 >> 1) / (threads*blocks)] elements of the current block
	uint64 id, i, j, minb, min, max, prime;
	uint64 global_min = completed * m + 1;


	id = threadIdx.x;
	min = global_min + (id * seglen);
	max = min + seglen - 2;

	//printf("Kernel %llu handles %11llu->%11llu\n", id, min, max);//works correctly

	prime = P[0];
	for (i = 1;(prime*prime <= max) && (i < plen);i++)
	{
		minb = ((min / prime) * prime);
		if (minb < min) minb += prime;
		if (~minb & 1)	minb += prime;
		for (j = minb;j <= max;j += (prime << 1))
			bitSet(mark, j - global_min + 1);
		prime = P[i];
	}

	//Print last found prime for last segment
	/*if (id == (blockDim.x - 1))
	{
		for (j = max; j >= min; j -= 2)
		{
			if (isPrime(mark, j - global_min + 1))
			{
				printf("Kernel %llu: %llu|%llu\n", id, j, max);
					break;
			}
		}
	}*/
}
////////////////////////////KERNELS END////////////////////////////

//     SEGMENTED SIEVE
//	n		RAM	    Time
// E07	   552KB   0.026s
// E08	   620KB   0.206s
// E09	   704KB   1.895s
// E10	   668KB   20.02s 
// E11     904KB   205.2s

//PARALLEL SEGMENTED SIEVE
//	n		RAM	    Time
// E10	      
// E11      

//Stats logged via Visual Studio Performance Profiler on i7 4790K @4.00GHz w/ 16GB DDR3 RAM and GTX 1070Ti

//Driver function
int main(uint32 argc, char* argv[])
{
	//Range: Data-type dependent
	uint64 n, m;
	printf("Enter n: ");
	scanf("%llu", &n);

	bound(n, m);
	uint32 threadCount = threadCalc(m);
	uint32 seglen = m/threadCount;

	uint32 plen = 0;
	uint32 *P = NULL;

	P = segboolsieve(n, m, plen);
	if (P == NULL)
	{
		printf("Memory Allocation Failure!\n");
		exit(0);
	}
	else
		printf("Last prime in utility sieve: %u @ index [%u]\n", P[plen - 1], plen - 1);

	uint32 segments = (uint32)((n + 1) / (m+1));		//No. of segments
	uint32 completed = 1;
	printf("\n%u segments(s) for [%llu->%llu]\n", segments - 1, m + 1, n);

	uint32 *dP;
	uint64 *dmark;
	uint64 *mark = (uint64 *)calloc((m / 128), sizeof(uint64));

	//Log execution time
	float GPU_TIME = 0.0;
	float temp_t;

	//CUDA Malloc
	cudaMalloc(&dP, (plen + 1) * (size));
	cudaMalloc(&dmark, ((m / 128) * size));

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	dim3 TPB(threadCount, 1, 1);

	cudaMemcpy(dP, P, plen * size, cudaMemcpyHostToDevice);

	while (completed < segments)
	{
		cudaMemcpy(dmark, mark, (m / 128) * size, cudaMemcpyHostToDevice);
		cudaEventRecord(start);
		SieveBlock << <1, TPB >> > (dP, dmark, m, completed, plen, seglen);
		cudaEventRecord(stop);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&temp_t, start, stop);
		//cudaMemcpy(dmark, mark, (m / 128) * size, cudaMemcpyHostToDevice);

		GPU_TIME += temp_t;
		completed++;
	}

	free(P);
	///free(mark);
	cudaFree(dP);
	///cudaFree(dmark);

	GPU_TIME /= 1000;
	printf("COMPUTE-PHASE GPU Time: %0.3f seconds\n", GPU_TIME);
	return 0;
}

