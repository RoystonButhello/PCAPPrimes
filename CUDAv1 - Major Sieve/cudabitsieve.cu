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
const uint64 pow2_32 = 4294967296;
const uint32 threads = 256;
const uint32 blocks = 8;
const uint64 segsize = pow2_32 / (threads*blocks);

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

__global__ void SieveBlock(uint32 *P, uint32 completed, uint32 plen)
{
	//Each thread sieves [(pow2_32 >> 1) / threads] elements of the current block
	uint64 mark[(segsize / 128) + 1];
	uint64 id, i, j, minb, min, max, prime, temp1, temp2;
	id = blockIdx.x*blockDim.x + threadIdx.x;
	min = (completed*pow2_32) + (id*segsize) + 1;
	max = min + segsize - 2;

	uint64 max_sieved = 0;
	bool max_set = false;
	for (i = 0;((uint64)P[i] * (uint64)P[i] <= max) && (i<plen);i++)
	{
		prime = P[i];
		minb = ((min / prime) * prime);
		if (minb < min) minb += prime;
		if (~minb & 1)	minb += prime;
		for (j = minb;j <= max;j += (prime << 1))
		{
			bitSet(mark, j - min + 1);
		}
		if (j - (prime << 1) > max_sieved) max_sieved = j - (prime << 1);
	}
	
	for (j = max; j >= min; j -= 2)
	{
		if (isPrime(mark, j - min + 1))
		{
			printf("Kernel %llu: %llu|%llu|%llu|%d\n", id , j, max_sieved, max, max_set);
			break;
		}
	}
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

	bool smallsieve = false;	//Use serial sieve for n<2^32
	if (n <= pow2_32)
	{
		smallsieve = true;
		printf("Rounded %llu to ", n);
		m = (uint64)sqrt(n);
		n = m * m;
		printf("%llu\n", n);
	}
	else if (n % pow2_32 > 0)	//If n>2^32 then round n to nearest multiple of 2^32
	{
		printf("Rounded %llu to ", n);
		n = ((n / pow2_32) + 1) * pow2_32;
		printf("%llu\n", n);
		m = (uint64)(sqrt(n));
	}

	uint32 plen = 0;
	uint32 *P = NULL;
	if (~n & 1) n--;
	if (~m & 1) m--;

	P = segboolsieve(n, m, plen, smallsieve);
	if (P == NULL)
	{
		printf("Memory Allocation Failure!\n");
		exit(0);
	}
	else
		printf("Last prime in utility sieve: %u @ index [%u]\n", P[plen - 1], plen - 1);

	if (smallsieve)
	{
		free(P);
		return 0;
	}

	uint32 chunkcount = (uint32)((n + 1) / pow2_32);		//No. of chunks
	uint32 completed = 1;
	printf("\n%u chunk(s) for [%llu->%llu]\n", chunkcount - 1, pow2_32 + 1, n);

	uint32 *dP;

	//Log execution time
	float GPU_TIME = 0.0;
	float temp_t;

	//CUDA Malloc
	cudaMalloc(&dP, (plen + 1) * (size));

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	dim3 TPB(threads, 1, 1);
	dim3 BPG(blocks, 1, 1);

	cudaMemcpy(dP, P, plen * size, cudaMemcpyHostToDevice);

	while (completed < chunkcount)
	{
		cudaEventRecord(start);
		SieveBlock << <BPG, TPB >> > (dP, completed, plen);
		cudaEventRecord(stop);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&temp_t, start, stop);

		GPU_TIME += temp_t;
		completed++;
	}

	free(P);
	cudaFree(dP);

	GPU_TIME /= 1000;
	printf("COMPUTE-PHASE GPU Time: %0.3f seconds\n", GPU_TIME);
	return 0;
}

