#include <iostream>
#include"Array.h"
#include<ctime>
using namespace std;

void quick(int *p, int len)
{
	for (int i = 0; i < len; i++)
	{
		cout << *(p + i) << " ";
	}
}
int quickSort(int *p, int *q)
{
	if (p == q || p - q>0)return 0;
	int pivot = *p;
	int *head = p;
	int *end = q;
	int flag = 0;
	while (p<q)
	{
		if (flag == 0)
		{
			while (*q > pivot)
			{
				q--;
				if (q == p)
				{
					*head = *p;
					*p = pivot;
					flag = 1;
					break;
				}
			}
		}
		if (flag == 0)
		{
			while (*p < pivot || *p == pivot)
			{
				p++;
				if (q == p)
				{
					*head = *p;
					*p = pivot;
					flag = 1;
					break;
				}
			}
		}
		if (flag == 0)
		{
			int temp;
			temp = *p;
			*p = *q;
			*q = temp;
		}
	}
	int a = quickSort(head, p - 1);
	int b = quickSort(p + 1, end);
	return 0;
}

int main()
{
	srand(time(0));
	long len, range;
	cout << "Enter the length:";
	cin >> len;
	cout << "Enter the range:";
	cin >> range;
	int *p = new int[len];
	//cout << p << endl; 

	//cout << "Enter " << len << " number:" << endl;
	for (int i = 0; i < len; i++)
	{
		//cin >> *(p + i);
		*(p + i) = rand() % range;
	}
	quick(p, len);
	int *q = p + len - 1;
	int a = quickSort(p, q);
	cout << endl ;
	quick(p, len);
	delete[]p;
	p = NULL;
	cout << "\n\n***********************************************\n";

	Array1d array1(1);
	array1.show();
	Array1d *p1 = new Array1d(1);
	p1->show();
	delete p1;
	

	return 0;
}
