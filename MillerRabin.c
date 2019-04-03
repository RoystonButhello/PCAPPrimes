//Input: A value n to be checked if prime or not
//Output: n if prime, closest prime lesser than n otherwise
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
typedef unsigned __int64 integer;
integer bitlen = (integer)8 * sizeof(integer);
const int pn = 9, P[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23 };

//Returns (a*b) % mod
integer mulMod(integer a, integer b, integer mod)
{
	a %= mod;
	b %= mod;
	integer res = 0;
	while (b)
	{
		if (b & 1)	//If b is odd
			res = (res + a) % mod;
		b >>= 1;
		a = ((integer)a << 1) % mod;
	}
	return res;
}

//Returns (base^exp) % mod
integer powerMod(integer base, integer exp, integer mod)
{
	integer res = 1;
	while (exp)
	{
		if (exp & 1) res = mulMod(res, base, mod);
		exp >>= 1;
		base = mulMod(base, base, mod);
	}
	return res;
}

_Bool isPrime(integer n)
{
	for (int i = 0; i < pn; i++)
		if (n % P[i] == 0)
			return n == P[i];	//Handle cases where n E P
	if (n < P[pn - 1]) return 0;
	integer s = 0, t = n - 1;
	while (~t & 1)		//While t is even and greater than 0
	{t >>= 1; s++;}		//Halve t, increment s
	for (int i = 0; i < pn; i++)
	{
		integer pt = powerMod(P[i], t, n);
		if (pt == 1) continue;
		_Bool flag = 0;
		for (int j = 0; j < s && !flag; j++)
		{
			if (pt == n - 1) flag = 1;
			pt = mulMod(pt, pt, n);
		}
		if (!flag) return 0;
	}
	return 1;
}

int main()
{
	//Measurements:
	//	  MILLER-RABIN ALGO
	//	n		RAM		  Time
	// 
	integer n;
	printf("Enter limit: "); scanf("%llu", &n);
	integer ceil = (integer)-1;
	ceil = ceil - 1 << (bitlen / 2 + 1) - 1;	//Avoid overflow
	if (n > ceil) n = ceil;
	if (~n & 1) n--;

	clock_t START_TIME, END_TIME;
	double  CPU_TIME;

	START_TIME = clock();
	if (!isPrime(n))
	{
		for (integer i = n ;i > 0;i-=2)
			if (isPrime(i)) { n = i; break; }
	}
	END_TIME = clock();

	printf("Result: %llu\n", n);
	CPU_TIME = ((double)(END_TIME - START_TIME)) / CLOCKS_PER_SEC;
	printf("\nCPU Time: %f seconds\n", CPU_TIME);
	return 0;
}
