#include  "Teacher.h"
#include <iostream>
#include <string>
using namespace std;

Teacher::Teacher(string _name, int _age)
{
	name = _name;
	age = _age;
	cout << "hi,i am in the constructor function" << endl;
}
Teacher::~Teacher()
{
	cout << "object has been destroyed" << endl;
}

void Teacher::setName(string _name)
{
	name = _name;
}
string Teacher::getName()
{
	return name;
}