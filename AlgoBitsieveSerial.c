//Algorithm 4.1 - Sequential Prime Sieve Using Bitwise Representation for Space Efficiency
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
//Key Data-type
typedef unsigned __int64 integer;	//Primary data-type; modifiable to allow flexibility
integer size = sizeof(integer);

//Utility function to calculate postive integer-powers
integer power(integer val, integer exp)
{
	integer temp = val;
	for (integer i = 1; i < exp; i++)
		temp *= val;
	if (temp == 0) temp--;	//Handle overflow when using types smaller than 64-bit
	return temp;
}

//Utility function to bound seed value to avoid overflow
integer bound(integer n)
{
	integer bitlen = size * 8;
	integer limit = power(2, bitlen) - power(2, (bitlen / 2) + 1) - 1;
	if (n > limit)
	{
		n = limit;
		printf("Bounded down to %llu\n", limit);
	}
	return n;
}

//Key-function to check if bit corresponding to a number <x> is set as composite or not
_Bool isPrime(integer *C, integer x)
{
	/*
	LOGIC via EXAMPLE
	Given a 32-bit integer, we have 2*size = 64 numbers covered per-byte of the array
	If x = 71, x/(2*size = 71/64 = 1 (71 represented by the the th bit of the 2nd element of C)
	x is halved (as even numbers are excluded) and undergoes bitwise AND with 31(0b11111) which is (size-1)
	In this case, (71/2) = 35
	35 AND 31 => 0b100101 AND 0b011111  = 0b101 => 5
	We now perform bitwise AND: C[1] AND (2^5)
	Essentially selecting the (5+1)th bit from C[1] and returning its value
	If 1, it is composite. Else it is prime
	C[0] should represent the numbers [2,3,5,7,9,...,63]
	C[1] should represent the numbers [65,67,69,...,127]
	*/
	return (C[x / (size << 4)] & ((integer)1 << ((x >> 1) & (size * 8 - 1)))) ? 0 : 1;
	/*
	integer index = x / (size << 4);
	integer byte = C[index];
	integer bitdex = (x >> 1) & (size * 8 - 1);
	integer bitval = ((integer)1 << bitdex);
	integer result = byte & bitval;
	return (result>0) ? 0 : 1;
	*///Debug block to track intermediate values
}

//Key-function to set bit corresponding to a number <x>
void mark(integer *C, integer x)
{
	//Set its corresponding bit in C[] to mark it as composite via bitwise OR
	C[x / (size << 4)] |= ((integer)1 << ((x >> 1) & (size * 8 - 1)));
	/*
	integer index = x / (size << 4);
	integer bitdex = (x >> 1) & (size * 8 - 1);
	integer orbyte = ((integer)1 << bitdex);
	integer mid = C[index] | orbyte;
	///printf("  %llu marked off at bitdex %llu of index %llu\n", x, bitdex,index); //Debug line
	C[index] = mid;
	*///Debug block to track intermediate values
}

//Driver function
int main()
{
	//Range: [Depends on defined data-type]

	//Statistics
	//	   BOOLEAN SIEVE
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

	//Stats logged via Visual Studio Performance Profiler on i7 4790K @4.00GHz w/ 16GB DDR3 RAM
	//Roughly 18% CPU usage for boolean array implementation and 13% CPU usage for bitwise sieve implementation
	//Can't take n>E10 due to ridiculous amount of time taken

	integer n;
	printf("Enter limit: ");
	scanf("%llu", &n);
	integer m = bound(n);
	//integer ceil = power(2, 32);
	//m = (n > ceil) ? ceil : n;				//Calculate primes upto m using phase 1, and primes>m using phase 2	
	integer len = (m / (2 * size * 8)) + 1;		//Only odd-values stored; 1 value per-bit of primary data-type
	if (m % 2 == 0) m--;

	//Boolean array initialized to false
	integer *C = (integer *)calloc(len, size);	//Represents [2,3,5,7,9,11,...] on per-bit basis
	if (C == NULL) { printf("Memory Allocation Failed!\n"); exit(1); }
	integer i, j;

	//Log execution time
	clock_t START_TIME, END_TIME;
	double  CPU_TIME;

	START_TIME = clock();
	for (i = 3; i*i <= m; i += 2)	//i->[3,5,7,9...,sqrt(n)] | i corresponds to C[(i-3)/2]
	{
		///printf("%llu selected at bitdex %llu of byte %llu\n", i, ((i >> 1) & (size * 8 - 1)), i / (size << 16),C[i / (size << 16)]); //Debug line
		if (isPrime(C, i))			//false implies <i> is prime
		{
			for (j = i * i; j <= m; j += (i << 1))	//j->[i^2,n] | increments by 2*i
			{
				mark(C, j);
				///printf(" %llu marked off by %llu\n", j, i); //Debug line
			}
			///if(i<100) printf("%llu stopped at %llu\n", i, j - 2 * i);	//Debug line
		}
	}
	END_TIME = clock();

	/*
	//Print all primes found for debugging
	short count = 1;
	printf("\nPrimes in range:\n2\t");
	for(i=3; i<n; i+=2)
	{
		if (isPrime(C,i))	  {printf("%llu\t",i); count++;}
		if (count%10 == 0) { printf("\n"); count /= 10;}
	}
	*/

	printf("\nLast 10 primes:\n");
	short count = 0;
	for (i = m; count < 10; i -= 2)
	{
		if (isPrime(C, i))
		{
			count++;
			printf("%llu\n", i);
		}
	}

	//Detail execution time
	CPU_TIME = ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;
	printf("\nCPU Time: %0.3f seconds\n", CPU_TIME);
	free(C);
	return 0;
}