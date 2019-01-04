#include <time.h>
#include <iostream>
#include <math.h>
#include <stdlib.h> 
#include <fstream>
#include <random>
#include <Continuous_Caster.h>
#include <Temperature.h>
#include <Steel.h>
#include <Gradientbasedalgorithm.h>
#include <Mesh.h>
#include <map> 
#include <string>  
#include <vector>
#include <Parameters_Control.h>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;
int main()
{
	const int section = 12, coolsection = 8, moldsection = 4;
	float ccml[section + 1] = { 0.0f,0.2f,0.4f,0.6f,0.8f,1.0925f,2.27f,4.29f,5.831f,9.6065f,13.6090f,19.87014f,28.599f };
	float h_init[section] = { 1380.0f,1170.0f,980.0f,800.0f,843.16f,605.05f,404.32f,342.83f,308.94f,281.64f,246.16f,160.96f };
	float h_upper[coolsection] = { 1343.16f, 855.05f, 654.32f, 492.83f, 428.94f, 401.64f, 386.16f, 260.96f };
	float h_lower[coolsection] = { 643.16f,455.05f, 254.32f, 242.83f, 228.94f, 181.64f, 146.16f, 100.96f };
	float taimmeantemperature[coolsection] = { 950.324f, 940.092f, 931.708f, 902.407f, 873.377f, 824.248f, 773.636f, 762.663f };
	float lx = 0.25f, ly = 28.599f, lz = 0.25f, tf = 4000.0f, vcast = -0.03f, rangeh = 50.0f, castingtemperature = 1558.0f, dh = 1.0f;
	float roll_radius = 0.095f;
	float roll_position[] = { 1.0925f, 1.3f, 1.52f, 1.74f, 1.99f, 2.27f, 2.55f, 2.84f, 3.13f, 3.42f, 3.71f, 4.0f, 4.29f, 4.595f, 4.904f, 5.213f, 5.522f, 5.831f,\
		6.143f, 6.460f, 6.773f, 7.07149f, 7.38449f, 7.70149f, 8.04849f, 8.36549f, 8.67849f, 8.97649f, 9.28949f, 9.60649f, 9.95199f, 10.28499f, \
		10.617f, 10.942f, 11.27599f, 11.60899f,  11.95199f, 12.28499f, 12.61799f, 12.94299f, 13.27599f, 13.60899f, 13.95199f, 14.28499f, 14.61799f,\
		14.94299f, 15.27599f, 15.60899f, 15.93393f, 16.23393f, 16.53393f, 16.83693f, 17.13993f, 17.44293f, 17.74593f, 18.05213f, 18.35513f, 18.65813f,\
		18.96113f, 19.26413f, 19.567f, 19.870f, 20.170f, 20.442f, 20.714f, 20.986f, 21.258f, 21.530f, 21.802f, 22.074f, 22.345f, 22.617f,\
		22.889f, 23.161f, 23.433f, 23.705f, 23.977f, 24.249f, 24.520f, 24.792f, 25.064f, 25.336f, 25.608f, 25.880f, 26.152f, 26.424f, \
		26.695f, 26.967f, 27.239f, 27.239f, 27.511f, 27.783f, 28.055f, 28.327f, 28.599f, 28.870f };
	int nx = 25, ny = 3001, nz = 25, tnpts = 20001;

	map<string, float> steel_componment;
	steel_componment["C"] = 0.12f;
	steel_componment["Mn"] = 0.35f;
	steel_componment["Si"] = 0.2f;
	steel_componment["S"] = 0.04f;
	steel_componment["P"] = 0.03f;
	steel_componment["Cr"] = 0.0f;
	steel_componment["Ni"] = 0.0f;
	steel_componment["Cu"] = 0.0f;
	steel_componment["Al"] = 0.0f;
	steel_componment["Ti"] = 0.0f;

	int predictstep = 2;
	Continuous_Caster CasterOne = Continuous_Caster(section, coolsection, moldsection, roll_radius, ccml, roll_position);
	Steel steel(steel_componment);
	Mesh continuous_casting_mesh = Mesh(nx, ny, nz, tnpts, lx, ly, lz, tf);
	Temperature SteelTemperature3dplant = Temperature(continuous_casting_mesh, vcast, CasterOne, steel);
	SteelTemperature3dplant.initcondition3d(castingtemperature);
	SteelTemperature3dplant.setcastingtemperature(castingtemperature);
	steel.print();
	
	Gradientbasedalgorithm Gradientmethod = Gradientbasedalgorithm(CasterOne, taimmeantemperature);

	int const n_regression = 3;
	const float vcast_n_regression[n_regression] = {-0.03f, -0.032f, -0.035f};
	vector<int> vcast_variation_tnpts = {0, 10001, 15001, 20001};

	Parameters_Control parameters_control = Parameters_Control(CasterOne, n_regression);
	parameters_control.set_vcast_mat(vcast_n_regression);
	clock_t t_start = clock();
	while (SteelTemperature3dplant.tstep < tnpts)
	{
		SteelTemperature3dplant.setvcast(n_regression, vcast_n_regression, vcast_variation_tnpts);
		if (SteelTemperature3dplant.tstep % 10 == 0)
			{
				Gradientmethod.computensparsejacobian(SteelTemperature3dplant, predictstep, h_init, dh);
				Gradientmethod.computegradient();
				Gradientmethod.computesearchdirection();
				Gradientmethod.linesearch();
				Gradientmethod.updateh(h_init, h_upper, h_lower);
				//Gradientmethod.print();
			}
			SteelTemperature3dplant.differencecalculation3d(h_init, predictstep);
			SteelTemperature3dplant.computemeantemperature3d();
			parameters_control.set_h_mat(h_init, vcast_variation_tnpts, SteelTemperature3dplant.tstep, Gradientmethod.costvalue);
			if (SteelTemperature3dplant.tstep % 50 == 0)
			{
				cout << "model timestep = " << SteelTemperature3dplant.tstep << " " << endl;
				SteelTemperature3dplant.print3d();
				Gradientmethod.outputdata(SteelTemperature3dplant, h_init);
			}
	}
	
	parameters_control.compute_parameters_values();
	parameters_control.print();
	parameters_control.output_parameters();

	clock_t t_end = clock();
	cout << "The running time is " << (t_end - t_start) << " (ms)" << endl;
}