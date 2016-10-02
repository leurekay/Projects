#include <iostream>
#include <string>
using namespace std;

class Teacher
{
public:
	Teacher(string name , int age );
	~Teacher();
	void setName(string _name);
	string getName();
	void setAge(int _age);
	int getAge();
private:
	string name;
	int age;
};