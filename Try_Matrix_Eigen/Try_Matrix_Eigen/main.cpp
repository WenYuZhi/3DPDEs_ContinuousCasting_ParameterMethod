#include <Eigen/Dense>
#include <iostream>
using namespace std;
using namespace Eigen;
int main()
{
	MatrixXd m = MatrixXd::Constant(3, 3, 1.0);
	MatrixXd n = MatrixXd::Constant(3, 3, 1.0);
	cout << m * n << endl;
	cout << m.rows() << endl;
	cout << m.cols() << endl;
}