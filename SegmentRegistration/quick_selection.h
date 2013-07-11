#include <iostream>
using namespace std;

template<class T>
int partition(T* input, int p, int r)
{
	T pivot = input[r];

	while ( p < r )
	{
		while ( input[p] < pivot )
			p++;

		while ( input[r] > pivot )
			r--;

		if ( input[p] == input[r] )
			p++;
		else if ( p < r ) {
			T tmp = input[p];
			input[p] = input[r];
			input[r] = tmp;
		}
	}

	return r;
}

template<class T>
T quick_select(T* input, int p, int r, int k)
{
	if ( p == r ) return input[p];
	if(k > (r-p+1)) return quick_select(input, p, r, r-p+1);
	int j = partition<T>(input, p, r);
	int length = j - p + 1;
	if ( length == k ) return input[j];
	else if ( k < length ) return quick_select<T>(input, p, j - 1, k);
	else  return quick_select<T>(input, j + 1, r, k - length);
}