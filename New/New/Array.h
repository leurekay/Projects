class Array1d
{
public:
	Array1d(int len);
	void show();
	int quickSort();

private:
	int m_iLen;
	int *m_pHead;
	int *m_pEnd;
};