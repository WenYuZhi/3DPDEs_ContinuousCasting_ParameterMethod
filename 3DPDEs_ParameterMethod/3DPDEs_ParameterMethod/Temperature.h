#ifndef TEMPERATURE_H
#define TEMPERATURE_H
#include <Steel.h>
#include <Continuous_Caster.h>
#include <Mesh.h>
#include <vector>
using namespace std;
class Temperature
{
    private:
	    float castingtemperature;
	    float h;
	    float* T_New;
	    float* T_Last;
	    bool disout;
	    
    public:
		const Continuous_Caster* mCasterOne;
		const Mesh* mesh;
		Steel* steel;
		float* T_Surface;
	    int tstep;
		float vcast;
	    float* meantemperature;
	    Temperature(const Mesh &, float, const Continuous_Caster &, Steel &);
	    Temperature(const Temperature &);
	    ~Temperature();
	    void differencecalculation3d(const float*, const int);
	    float boundarycondition3d(const float*, int);
		float boundarycondition_with_roll(const float*, int);
	    void initcondition3d(float);
	    void initcondition3d(float*);
	    void print3d(int);
	    void print3d();
	    void computemeantemperature3d();
	    void setvcast(const int, const float*, vector<int>);
		void setcastingtemperature(float);
	    void operator=(const Temperature &);
	    friend class Gradientbasedalgorithm;
};
#endif