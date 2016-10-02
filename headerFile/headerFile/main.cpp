#include <iostream>
#include <string>
#include  "Teacher.h"
using namespace std;



int main()
{
	Teacher t1("Allen", 17);
	//t1.setName("Bill");
	cout << t1.getName()<<endl;
	return 0;
}