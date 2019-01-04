#include <Continuous_Caster.h>
#include <iostream>
using namespace std; 
Continuous_Caster::Continuous_Caster(int section, int coolsection, int moldsection, float roll_radius, float* ccml, float roll_position [])
{
	this->section = section;
	this->coolsection = coolsection;
	this->moldsection = moldsection;
	this->length_roll_contact = 0.1f * roll_radius;
	this->roll_radius = roll_radius;
	this->n_roll = sizeof(roll_position) / sizeof(float);
	this->ccml = new float[section + 1];
	for (int i = 0; i < section + 1; i++)
		this->ccml[i] = ccml[i];
	this->roll_position = new float[this->n_roll];
	for (int i = 0; i < this->n_roll; i++)
		this->roll_position[i] = roll_position[i];
}

Continuous_Caster::~Continuous_Caster()
{
	delete[] ccml;
	delete[] roll_position;
}

void Continuous_Caster::print()
{
	cout << "section = " << section << " ";
	cout << "coolsection = " << coolsection << " ";
	cout << "moldsection = " << moldsection << " " << endl;
	for (int i = 0; i < section; i++)
		cout << ccml[i] << ", ";
}