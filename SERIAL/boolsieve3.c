//Algorithm 4.1 - Sequential Prime Sieve Using Fast Marking w/ Bitwise Representation for Space Efficiency
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
	for (integer i = 1; i < exp; i++) temp *= val;
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
	return (C[x / (size<<1)] & (1 << ((x >> 1) & (size - 1)))) ? 0 : 1;
}

//Key-function to set bit corresponding to a number <x>
void mark(integer *C, integer x)
{
	//Set its corresponding bit in C[] to mark it as composite via bitwise OR
	C[x / (size << 1)] |= (1 << ((x >> 1) & (size - 1)));
	//printf("%llu marked off\n", x); //Debug line
}
//Driver function
int main()
{
	//Range: [Depends on defined data-type]
	//Space efficiency: When indexing using 32-bit numbers, it uses about 2GB of RAM to manage numbers upto 4,294,967,296
	//					When indexing using 64-bit numbers, it uses about 4GB of RAM to manage numbers upto 86,580,000,000
	//86580000000
	integer n;
	printf("Enter limit: ");
	scanf("%llu", &n);
	n = bound(n);
	integer len = n / (2 * size);		  //Only odd-values stored; 1 value per-bit of primary data-type
	if (n % 2 == 0) n--;

	//Boolean array initialized to false
	integer *C = (integer *)calloc(len, size);	//Represents [2,3,5,7,9,11,...] on per-bit basis
	if (C == NULL) { printf("Memory Allocation Failed!\n"); exit(1); }
	integer i, j;

	//Log execution time
	clock_t START_TIME, END_TIME;
	double  CPU_TIME;

	START_TIME = clock();
	for (i = 3; i*i <= n; i += 2)	//i->[3,5,7,9...,sqrt(n)] | i corresponds to C[(i-3)/2]
	{
		if (isPrime(C, i))			//false implies <i> is prime
		{
			for (j = i * i; j <= n; j += (i<<1))	//j->[i^2,n] | increments by 2*i
				mark(C, j);
			//printf("%llu stopped at %llu\n", i, j - 2 * i);
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
	for (i = n; count < 10; i -= 2)
	{
		if (isPrime(C, i))
		{
			count++;
			printf("%llu\n", i);
		}
	}
	
	
	integer pcount = 0;
	for (i = 0, pcount = 0; i < len; i++)
		if (isPrime(C,i)) pcount++;
	printf("%llu primes found\n", pcount);
	
	//Detail execution time
	CPU_TIME = ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;
	printf("\nCPU Time: %0.3f seconds\n", CPU_TIME);
	return 0;
}