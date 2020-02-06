#include <vector>
#include <iostream>
#include "updateable_heap.h"
using namespace std;

ostream& operator<<(ostream &str, heapElemT he)
{
	str << *(he.left);
	return str;
}


int operator>(const heapElemT &a, const heapElemT &b) // compares the probability log ratios
{
	return ( (*a.left)->probLogRatio > (*b.left)->probLogRatio );
}
