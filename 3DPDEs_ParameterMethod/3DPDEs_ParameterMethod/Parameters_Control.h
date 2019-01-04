#ifndef PARAMETERS_CONTROL_H
#define PARAMETERS_CONTROL_H
#include <Temperature.h>
#include <Continuous_Caster.h>
#include <Eigen/Dense>
#include <vector>
using namespace std;
using namespace Eigen;
class Parameters_Control
{
    public:
		MatrixXd vcast;
		MatrixXd h_static;
		MatrixXd parameters_control_coeff;
		int coolsection, moldsection, n_dim;
	    Parameters_Control(const Continuous_Caster &, int);
	    void compute_parameters_values();
		void set_vcast_mat(const float*);
		void set_h_mat(const float*, vector<int>, int, float);
		void output_parameters();
		void print();
};
#endif
