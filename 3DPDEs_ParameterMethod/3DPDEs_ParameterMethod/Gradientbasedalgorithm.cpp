#include <Gradientbasedalgorithm.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h> 
#include <time.h>  
using namespace std;
Gradientbasedalgorithm::Gradientbasedalgorithm(Continuous_Caster & CasterOne, float* taimmeantemperature)
{
	mCasterOne = &CasterOne;
	coolsection = mCasterOne->coolsection;
	costvalue = 0.0f;
	allmeantemperature = new float*[mCasterOne->coolsection];
	for (int i = 0; i < mCasterOne->coolsection; i++)
		allmeantemperature[i] = new float[coolsection];

	for (int i = 0; i < coolsection; i++)
		for (int j = 0; j < coolsection; j++)
			allmeantemperature[i][j] = 0.0f;

	Jacobian = new float*[mCasterOne->coolsection];
	for (int i = 0; i < mCasterOne->coolsection; i++)
		Jacobian[i] = new float[coolsection];

	for (int i = 0; i < coolsection; i++)
		for (int j = 0; j < coolsection; j++)
			Jacobian[i][j] = 0.0f;

	staticmeantemperature = new float[coolsection];
	for (int i = 0; i < coolsection; i++)
		staticmeantemperature[i] = 0.0f;

	gradient = new float[coolsection];
	for (int i = 0; i < coolsection; i++)
		gradient[i] = 0.0f;

	dk = new float[coolsection];
	for (int i = 0; i < coolsection; i++)
		dk[i] = 0.0f;

	this->taimmeantemperature = new float[coolsection];
	for (int i = 0; i < mCasterOne->coolsection; i++)
		this->taimmeantemperature[i] = taimmeantemperature[i];
}

Gradientbasedalgorithm::~Gradientbasedalgorithm()
{
	for (int i = 0; i < mCasterOne->coolsection; i++)
		delete[] allmeantemperature[i];
	delete[] allmeantemperature;
	for (int i = 0; i < mCasterOne->coolsection; i++)
		delete[] Jacobian[i];
	delete[] Jacobian;
	delete[] taimmeantemperature;
	delete[] staticmeantemperature;
	delete[] gradient;
	delete[] dk;
}


void Gradientbasedalgorithm::computejacobian(const Temperature & SteelTemperature3dmodel, int predictstep, float* h_init, float dh)
{
	Temperature SteelTemperature3dtemp = SteelTemperature3dmodel;
	float *h_perturbation = new float[mCasterOne->section];
	for (int j = 0; j < mCasterOne->section; j++)
		h_perturbation[j] = h_init[j];
	for (int i = 0; i < coolsection + 1; i++)
	{
		SteelTemperature3dtemp = SteelTemperature3dmodel;
		SteelTemperature3dtemp.computemeantemperature3d();

		if (i == coolsection)
		{
			SteelTemperature3dtemp.differencecalculation3d(h_init, predictstep);
			SteelTemperature3dtemp.computemeantemperature3d();
			for (int j = 0; j < coolsection; j++)
				staticmeantemperature[j] = SteelTemperature3dtemp.meantemperature[j];
		}
		else
		{
			for (int j = 0; j < mCasterOne->section; j++)
				h_perturbation[j] = h_init[j];
			h_perturbation[mCasterOne->moldsection + i] += dh;
			SteelTemperature3dtemp.differencecalculation3d(h_perturbation, predictstep);
			SteelTemperature3dtemp.computemeantemperature3d();
			for (int j = 0; j < coolsection; j++)
				allmeantemperature[i][j] = SteelTemperature3dtemp.meantemperature[j];
		}
	}

	for (int i = 0; i < coolsection; i++)
		for (int j = 0; j < coolsection; j++)
			Jacobian[i][j] = (allmeantemperature[i][j] - staticmeantemperature[j]) / dh;
}

void Gradientbasedalgorithm::computesparsejacobian(const Temperature & SteelTemperature3dmodel, int predictstep, float* h_init, float dh)
{
	Temperature SteelTemperature3dtemp = SteelTemperature3dmodel;
	float *h_perturbation = new float[mCasterOne->section];
	for (int j = 0; j < mCasterOne->section; j++)
		h_perturbation[j] = h_init[j];
	int n_sparse = 1;

	SteelTemperature3dtemp.differencecalculation3d(h_init, predictstep);
    SteelTemperature3dtemp.computemeantemperature3d();

	for (int j = 0; j < coolsection; j++)
		staticmeantemperature[j] = SteelTemperature3dtemp.meantemperature[j];

    for (int j = 0; j < mCasterOne->section; j++)
		h_perturbation[j] = h_init[j];
	for (int j = 0; j < mCasterOne->coolsection; j++)
	    h_perturbation[mCasterOne->moldsection + j] += dh;

	SteelTemperature3dtemp = SteelTemperature3dmodel;
	SteelTemperature3dtemp.differencecalculation3d(h_perturbation, predictstep);
	SteelTemperature3dtemp.computemeantemperature3d();
	for (int m = 0; m < coolsection; m++)
		for (int j = 0; j < coolsection; j++)
			if (m == j)
				 allmeantemperature[m][j] = SteelTemperature3dtemp.meantemperature[j];

	for (int i = 0; i < coolsection; i++)
		for (int j = 0; j < coolsection; j++)
			if (i == j)
				Jacobian[i][j] = (allmeantemperature[i][j] - staticmeantemperature[j]) / dh;
			else
				Jacobian[i][j] = 0.0f;
}

void Gradientbasedalgorithm::computensparsejacobian(const Temperature & SteelTemperature3dmodel, int predictstep, float* h_init, float dh)
{
	Temperature SteelTemperature3dtemp = SteelTemperature3dmodel;
	float *h_perturbation = new float[mCasterOne->section];
	for (int j = 0; j < mCasterOne->section; j++)
		h_perturbation[j] = h_init[j];
	int n_sparse = 2;

	SteelTemperature3dtemp.differencecalculation3d(h_init, predictstep);
	SteelTemperature3dtemp.computemeantemperature3d();
	for (int j = 0; j < coolsection; j++)
		staticmeantemperature[j] = SteelTemperature3dtemp.meantemperature[j];

	for (int j = 0; j < mCasterOne->section; j++)
		h_perturbation[j] = h_init[j];
	for (int j = 0; j < coolsection; j++)
		if (j % n_sparse == 0)
			h_perturbation[mCasterOne->moldsection + j] += dh;

	SteelTemperature3dtemp = SteelTemperature3dmodel;
	SteelTemperature3dtemp.differencecalculation3d(h_perturbation, predictstep);
	SteelTemperature3dtemp.computemeantemperature3d();

	for (int j = 0; j < coolsection; j++)
		allmeantemperature[0][j] = SteelTemperature3dtemp.meantemperature[j];

	for (int j = 0; j < mCasterOne->section; j++)
		h_perturbation[j] = h_init[j];
	for (int j = 0; j < coolsection; j++)
		if (j % n_sparse == 1)
			h_perturbation[mCasterOne->moldsection + j] += dh;

	SteelTemperature3dtemp = SteelTemperature3dmodel;
	SteelTemperature3dtemp.differencecalculation3d(h_perturbation, predictstep);
	SteelTemperature3dtemp.computemeantemperature3d();

	for (int j = 0; j < coolsection; j++)
		allmeantemperature[1][j] = SteelTemperature3dtemp.meantemperature[j];

	for (int i = 0; i < coolsection; i++)
		for (int j = 0; j < coolsection; j++)
			if (i % n_sparse == 0 && (i == j || i == (j - 1)))
				Jacobian[i][j] = (allmeantemperature[0][j] - staticmeantemperature[j]) / dh;
			else if (i % n_sparse == 1 && (i == j || i == (j - 1)))
				Jacobian[i][j] = (allmeantemperature[1][j] - staticmeantemperature[j]) / dh;
			else
				Jacobian[i][j] = 0.0f;
	delete[] h_perturbation;
}

void Gradientbasedalgorithm::computestochasticnsparsejacobian(const Temperature & SteelTemperature3dmodel, int predictstep, float* h_init, float dh, int random_index)
{
	Temperature SteelTemperature3dtemp = SteelTemperature3dmodel;
	float *h_perturbation = new float[mCasterOne->section];
	for (int j = 0; j < mCasterOne->section; j++)
		h_perturbation[j] = h_init[j];
	int n_sparse = 2;

	SteelTemperature3dtemp.differencecalculation3d(h_init, predictstep);
	SteelTemperature3dtemp.computemeantemperature3d();
	for (int j = 0; j < coolsection; j++)
		staticmeantemperature[j] = SteelTemperature3dtemp.meantemperature[j];

	if (random_index == 0)
	{
		for (int j = 0; j < mCasterOne->section; j++)
			h_perturbation[j] = h_init[j];
		for (int j = 0; j < coolsection; j++)
			if (j % n_sparse == 0)
				h_perturbation[mCasterOne->moldsection + j] += dh;

		SteelTemperature3dtemp = SteelTemperature3dmodel;
		SteelTemperature3dtemp.differencecalculation3d(h_perturbation, predictstep);
		SteelTemperature3dtemp.computemeantemperature3d();

		for (int j = 0; j < coolsection; j++)
			allmeantemperature[0][j] = SteelTemperature3dtemp.meantemperature[j];

		for (int i = 0; i < coolsection; i++)
			for (int j = 0; j < coolsection; j++)
				if (i % n_sparse == 0 && (i == j || i == (j - 1)))
					Jacobian[i][j] = (allmeantemperature[0][j] - staticmeantemperature[j]) / dh;
	}
	
	else if (random_index == 1)
	{
		for (int j = 0; j < mCasterOne->section; j++)
			h_perturbation[j] = h_init[j];
		for (int j = 0; j < coolsection; j++)
			if (j % n_sparse == 1)
				h_perturbation[mCasterOne->moldsection + j] += dh;

		SteelTemperature3dtemp = SteelTemperature3dmodel;
		SteelTemperature3dtemp.differencecalculation3d(h_perturbation, predictstep);
		SteelTemperature3dtemp.computemeantemperature3d();

		for (int j = 0; j < coolsection; j++)
			allmeantemperature[1][j] = SteelTemperature3dtemp.meantemperature[j];

		for (int i = 0; i < coolsection; i++)
			for (int j = 0; j < coolsection; j++)
				if (i % n_sparse == 1 && (i == j || i == (j - 1)))
					Jacobian[i][j] = (allmeantemperature[1][j] - staticmeantemperature[j]) / dh;
	}

	else
		cout << "random index variable is not right";
	
	delete[] h_perturbation;
}

int Gradientbasedalgorithm::generaterandomindex()
{
	srand((unsigned)time(NULL));
	return int(rand() % 2);
}

void Gradientbasedalgorithm::computegradient()
{
	for (int i = 0; i < coolsection; i++)
	{
		gradient[i] = 0.0;
		for (int j = 0; j < coolsection; j++)
			gradient[i] += (taimmeantemperature[j] - staticmeantemperature[j]) * Jacobian[i][j];
	}
}

void Gradientbasedalgorithm::computesearchdirection()
{
	for (int i = 0; i < coolsection; i++)
	    dk[i] = gradient[i];
}

void::Gradientbasedalgorithm::linesearch()
{
	float step1 = 0.0f, step2 = 0.0f, eps = 0.01f;
	for (int i = 0; i < coolsection; i++)
		for (int j = 0; j < coolsection; j++)
		{
			step1 += (staticmeantemperature[i] - taimmeantemperature[i]) * Jacobian[i][j];
			step2 += Jacobian[i][j] * dk[j] * Jacobian[i][j] * dk[j];
		}
	if (step2 < eps)
		step = fabs(step1 / (step2 + eps));
	else
		step = fabs(step1 / step2);
}

void::Gradientbasedalgorithm::updateh(float*h_init, float*h_upper, float*h_lower)
{
	for (int i = 0; i < coolsection; i++)
	{
		h_init[i + mCasterOne->moldsection] = h_init[i + mCasterOne->moldsection] + step * dk[i];
		if (h_init[i + mCasterOne->moldsection] > h_upper[i])
			h_init[i + mCasterOne->moldsection] = h_upper[i];
		else if (h_init[i + mCasterOne->moldsection] < h_lower[i])
			h_init[i + mCasterOne->moldsection] = h_lower[i];
	}
	costvalue = 0.0;
	for (int i = 0; i < coolsection; i++)
		costvalue += (staticmeantemperature[i] - taimmeantemperature[i]) * (staticmeantemperature[i] - taimmeantemperature[i]);
	costvalue = float(pow(costvalue / coolsection, 0.5));
}


void Gradientbasedalgorithm::init(float**allmeantemperature, float *staticmeantemperature, float dh)
{
	this->dh = dh;
	for (int i = 0; i < coolsection; i++)
		for (int j = 0; j < coolsection; j++)
			this->allmeantemperature[i][j] = allmeantemperature[i][j];
	for (int i = 0; i < coolsection; i++)
		this->staticmeantemperature[i] = staticmeantemperature[i];
}

void::Gradientbasedalgorithm::validation(float*m_measuredtemperaturetemp, float*hinit, float*htemp)
{
	float msetemperature = 0.0;
	for (int i = 0; i < coolsection; i++)
		msetemperature += (m_measuredtemperaturetemp[i] - staticmeantemperature[i]) * (m_measuredtemperaturetemp[i] - staticmeantemperature[i]);

	float mseh = 0.0;
	for (int i = 0; i < coolsection; i++)
		mseh += (hinit[i + mCasterOne->moldsection] - htemp[i + mCasterOne->moldsection]) * (hinit[i + mCasterOne->moldsection] - htemp[i + mCasterOne->moldsection]);
	cout << "mseh = " << mseh << ", ";
	cout << "msetemperature = " << msetemperature << endl;
}

void Gradientbasedalgorithm::print()
{
	cout << endl;
	cout << "Jacobian = " << endl;
	for (int i = 0; i < coolsection; i++)
	{
	    for (int j = 0; j < coolsection; j++)
	        cout << Jacobian[i][j] << ",";
	    cout << endl;
	}

	cout << "staticmeantemperature = " << endl;
	for (int i = 0; i < coolsection; i++)
		cout << staticmeantemperature[i] << ", ";
	cout << endl;

	cout << "staticmeantemperature - taimmeantemperature = " << endl;
	for (int i = 0; i < coolsection; i++)
		cout << fabs(staticmeantemperature[i] - taimmeantemperature[i]) << ", ";
	cout << endl;

	cout << "Gradient = ";
	for (int i = 0; i < coolsection; i++)
		cout << gradient[i] << ", ";
	cout << endl;
	cout << "step = " << step << endl;
	cout << "costvalue = " << costvalue << endl;
}

void Gradientbasedalgorithm::outputdata(Temperature & m_SteelTemperature, float* h_init)
{
	ofstream outputfile;
	outputfile.open("C:\\T_Result_3DCPU_MPC.csv", ios::app);
	for (int i = 0; i < m_SteelTemperature.mesh->ny; i++)
		outputfile << m_SteelTemperature.T_Surface[i] << ",";
	outputfile << endl;
	outputfile.close();

	outputfile.open("C:\\Temperature3DCPU_MPC.csv", ios::app);
	outputfile << m_SteelTemperature.tstep * m_SteelTemperature.mesh->tao << "," << costvalue << ",";
	outputfile << m_SteelTemperature.vcast << ",";
	for (int i = 0; i < coolsection; i++)
		outputfile << h_init[i + mCasterOne->moldsection] << ",";
	for (int i = 0; i < coolsection; i++)
		outputfile << m_SteelTemperature.meantemperature[i] << ",";
	for (int i = 0; i < coolsection; i++)
		outputfile << (m_SteelTemperature.meantemperature[i] - taimmeantemperature[i]) << ",";
	for (int i = 0; i < coolsection; i++)
		outputfile << gradient[i] << ",";
	outputfile << endl;
	outputfile.close();
}