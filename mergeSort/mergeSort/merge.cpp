#include <iostream>
#include <stdlib.h>
using namespace std;
void bubbleSort(int len,int a[])
{
	for (int i = 0; i < len - 1; i++)
	{
		for (int j = 0; j < len - 1 - i; j++)
		{
			if (a[j] > a[j + 1])
			{
				int temp = a[j];
				a[j] = a[j + 1];
				a[j + 1] = temp;
			}
		}
	}
}

int main()
{	
	int a[] = {4,6,8,7,9,7,6,7,84,34,84,34,4,4,4,5,6,4,34,3,7,8,9,9,98,5,8,6,13,645,4,45,4,54,4,45,8,35,4,54,648,4454,452,45,4,5,
	52,54,24,2,54,42,42,452,24,40};
	int len = sizeof(a) / 4;
	bubbleSort(len,a);
	for (int i = 0; i < len; i++)cout << a[i]<<' ' ;
	cout << endl;
	cout <<"number="<<sizeof(a)/4<< endl;
	system("pause");
	return 0;
}