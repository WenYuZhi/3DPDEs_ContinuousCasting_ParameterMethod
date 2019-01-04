#ifndef CONTINUOUS_CASTER_H
#define CONTINUOUS_CASTER_H
class Continuous_Caster
{
    public:
	    int section, coolsection, moldsection, n_roll;
	    float *ccml, *roll_position;
		float roll_radius, length_roll_contact;
	    Continuous_Caster(int, int, int, float, float*, float []);
	    ~Continuous_Caster();
	    void print();
};
#endif

