#include <Parameters_Control.h>
#include <Eigen/Dense>
#include <Continuous_Caster.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;
using namespace Eigen;
Parameters_Control::Parameters_Control(const Continuous_Caster & CasterOne, int n_dim)
{
	this->coolsection = CasterOne.coolsection;
	this->moldsection = CasterOne.moldsection;
	this->n_dim = n_dim;
	vcast = MatrixXd::Constant(n_dim, 3, 1.0f);
	h_static = MatrixXd::Constant(n_dim, this->coolsection, 1.0f);
	parameters_control_coeff = MatrixXd::Constant(3, this->coolsection, 0.0f);
}

void Parameters_Control::set_vcast_mat(const float* vcast_n_regression)
{
	for (int i = 0; i < n_dim; i++)
	{
		vcast(i, 0) = (vcast_n_regression[i] * 60.0f) * (vcast_n_regression[i] * 60.0f);
		vcast(i, 1) = fabs(vcast_n_regression[i]) * 60.0f;
	}
}

void Parameters_Control::set_h_mat(const float* h_init, vector<int> vcast_variation_tnpts, int tstep, float costvalue)
{
	int static num = 0;
	vector<int>::iterator iter = find(vcast_variation_tnpts.begin(), vcast_variation_tnpts.end(), tstep + 1); //返回的是一个迭代器指针
	if (iter != vcast_variation_tnpts.end() && costvalue <= 11110.5f)
	{
		for (int i = 0; i < h_static.cols(); i++)
			h_static(num, i) = h_init[i + moldsection];
		num++;
	}
}

void Parameters_Control::compute_parameters_values()
{
	MatrixXd vcast_inv = MatrixXd::Constant(vcast.rows(), vcast.cols(), 0.0f);
	vcast_inv = vcast.inverse();
	cout << "vcast_inv = " << vcast_inv << endl;
	parameters_control_coeff = vcast_inv * h_static;
}

void Parameters_Control::output_parameters()
{
	ofstream outputfile;
	outputfile.open("C:\\parameters_abc.csv", ios::app);
	outputfile << parameters_control_coeff;
	outputfile.close();
}

void Parameters_Control::print()
{
	cout << "vcast_mat = " << vcast << endl;
	cout << "h_mat = " << h_static << endl;
	cout << "abc = " << parameters_control_coeff << endl;
}