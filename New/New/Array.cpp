#include<iostream>
#include "Array.h"
using namespace std;

Array1d::Array1d(int len)
{
	m_iLen = len;
	m_pHead = new int[m_iLen];
	m_pEnd = m_pHead + m_iLen - 1;
	
	cout << "\nenter " << m_iLen << " " << "number:\n";
	for (int i = 0; i < m_iLen; i++)
	{
		cin >> *(m_pHead + i);
	}
}

void Array1d::show()
{
	for (int i = 0; i < m_iLen; i++)
	{
		cout << *(m_pHead + i) << " ";
	}
}
/*
int Array1d::quickSort()
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
*/