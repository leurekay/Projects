#include <iostream>
using namespace std;
int main()
{
	const int months_of_year = 12;
	double principal = 0;
	cout << "enter the principal:";
	cin >> principal;
	double interest = 0.0;
	cout << "enter the interest ratio:";
	cin >> interest;
	int years = 0;
	cout << "int the years of loan:";
	cin >> years;
	double month_interest = interest / months_of_year;
	long months = years*months_of_year;
	cout << "principal:" << principal << " " << "month_interest:" << month_interest << endl;
	return 0;

}