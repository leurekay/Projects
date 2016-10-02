#include <iostream>
#include <string>
using namespace std;

class Coor
{
public:
	int x;
	int y;
	void printx()
	{
		cout<< x << endl;
	}
	void printy()
	{
		cout << y << endl;
	}
};

int main()
{
	Coor stackCoor;
	stackCoor.x = 10;
	stackCoor.y = 20;
	stackCoor.printx();

	Coor *p = new  Coor();
	p->x = 100;
	p->y= 200;
	p->printy();
	delete p;
	p = NULL;

	char a[4] = { 'a','b','c','d'};
	string hobby ("football");
	cout << a << endl;
	cout<<hobby<< endl;
	return 0;
}