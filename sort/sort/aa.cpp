#include<iostream>
using namespace std;   

void swap(int &x, int &y)
{
	int temp;
	temp = x;
	x = y;
	y = temp;
}

int main()
{	
	int i, j;
	cin >> i >> j;
	cout << i << "   "<<j << endl;
	swap(i, j);
	cout << i << "   " << j << endl;
	int a[][2]= { 0,1,2,3,4,5};
 
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 2; j++)
		{
			cout << "a[" << i << "]" << "[" << j << "]=" << a[i][j] << endl;
		}
	}
	cout << &i << endl;
	system("pause");
	return 0;
}