#include <iostream>
#include <string>
using namespace std;
class Student
{
public:
	void setName(string _name)
	{
		name = _name;
	}
	string getName()
	{ 
		return name; 
	}
	void setAge(int _age)
	{
		age = _age;
	}
	int getAge()
	{
		return age;
	}
private:
	string name;
	int score;
	int age;
};
int main()
{
	Student stu;
	stu.setName("Bill");
	stu.setAge(500);
	
	cout <<stu.getName()<<stu.getAge()<< endl;
	return 0;
}