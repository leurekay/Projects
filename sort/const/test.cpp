#include<iostream>
using namespace std;

int main()
{
	int x = 3;
	int y = 5;
	const int *p1 = &x;  //��
	int const *p2 = null;  //��
	p1 = &y;
	*p1 = 8;
	return 0;
}