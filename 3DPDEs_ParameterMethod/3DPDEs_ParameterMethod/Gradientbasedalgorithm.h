#ifndef GRADIENTBASEALGORITHM_H
#define GRADIENTBASEALGORITHM_H
#include <Continuous_Caster.h>
#include <Temperature.h>

class Gradientbasedalgorithm
{
    public:
	int coolsection;
	float** allmeantemperature;
	float* staticmeantemperature;
	float* taimmeantemperature;
	float** Jacobian;
	float* gradient;
	float* dk;
	float costvalue;
	float dh;
	float step;
	Continuous_Caster* mCasterOne;
	Gradientbasedalgorithm(Continuous_Caster &, float*);
	~Gradientbasedalgorithm();
	void computejacobian(const Temperature &, int, float*, float);
	void computesparsejacobian(const Temperature &, int, float*, float);
	void computensparsejacobian(const Temperature &, int, float*, float);
	void computestochasticnsparsejacobian(const Temperature &, int, float*, float, int);
	int generaterandomindex();
	void computegradient();
	void computesearchdirection();
	void init(float**, float*, float);
	void linesearch();
	void updateh(float*, float*, float*);
	void validation(float*, float*, float*);
	void print();
	void outputdata(Temperature &, float*);
};
#endif
