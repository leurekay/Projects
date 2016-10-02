#include <iostream>
using namespace std;

void quick(int *p,int len)
{
	for (int i = 0; i < len; i++)
	{
		cout << *(p+i)<<" ";
	}
}
int quickSort(int *p, int *q)
{	
	if (p == q)return 0;
	int pivot = *p;
	int *head = p;
	int *end = q;
	while (p<q) 
	{
		while (*q > pivot)
		{
			q--;
			if (q == p)
			{
				*head = *p;
				*p = pivot;
				break;
			}
		}
		while (*p < pivot || *p == pivot)
		{
			if (q == p)
			{
				*head = *p;
				*p = pivot;
				break;
			}
			p++;
			
		}
		int temp;
		temp = *p;
		*p = *q;
		*q = temp;
	}
	int a=quickSort(head, p - 1);
	int b=quickSort(p + 1, end);
	return 0;
}

int main()
{

	long len;
	cout << "Enter the length:";
	cin >> len;
	int *p= new int[len];
	cout << p << endl;

	cout << "Enter " << len << " number:" << endl;
	for (int i = 0; i < len; i++)
	{
		cin >> *(p+i);
	}
	quick(p, len);
	int *q = p + len - 1;
	int c=quickSort(p, q);
	cout << endl;
	quick(p, len);
	delete []p;
	p = NULL;
	return 0;
}