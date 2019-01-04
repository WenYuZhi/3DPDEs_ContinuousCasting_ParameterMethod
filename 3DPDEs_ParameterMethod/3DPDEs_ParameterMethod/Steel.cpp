#include <Steel.h>

Steel::Steel(map<string, float> steel_componment)
{
	this->steel_componment = steel_componment;
	pho = 7250.0f;
	ce = 660.0f;
	L = 265600.0f;
	lamda_liquid = 35.0f - 0.3574f * steel_componment["Cr"] - 0.5116f * steel_componment["Mn"] + 0.0014f * steel_componment["Ni"];
	lamda_solid = 20.76f - 3.2627f * steel_componment["C"] - 0.7598f * steel_componment["Si"] - 0.1432f * steel_componment["Mn"];
	temperature_liquid = 1536.6f - (88.0f * steel_componment["C"] + 8.0f * steel_componment["Si"] + 5.0f * steel_componment["Mn"] + 25.0f * steel_componment["Cu"] + 1.5f * steel_componment["Cr"] + 4.0f * steel_componment["Ni"] + 2.0f * steel_componment["Mn"] + 18.0f * steel_componment["Ti"]);
	temperature_solid = 1527.0f - (187.5f * steel_componment["C"] + 700.0f * steel_componment["S"] + 500.0f * steel_componment["P"] + 20.5f * steel_componment["Si"] + 11.5f*steel_componment["Ni"] + 6.5f * steel_componment["Mn"] + 5.5f * steel_componment["Al"] + 2.0f *steel_componment["Cr"]);
}

void Steel::physical_parameters(float T)
{
	float fs, enhance_factor = 2.0f;
	if (T < temperature_liquid)
	{
		fs = 0.0f;
		lamda = lamda_solid;
		ce = 660.0f;
	}

	else if (T > temperature_liquid)
	{
		fs = 1.0f;
		lamda = lamda_solid + 0.009f * T;
		ce = 660.0f;
	}

	else
	{
		fs = (T - temperature_solid) / (temperature_liquid - temperature_solid);
		lamda = fs * lamda_solid + enhance_factor * (1 - fs) * lamda_liquid;
		ce = ce + L / (temperature_liquid - temperature_solid);
	}
}

void Steel::print()
{
	map<string, float>::iterator iter;
	for (iter = steel_componment.begin(); iter != steel_componment.end(); iter++)
		cout << iter->first << ":" << iter->second << "%  ";
	cout << endl;
	cout << "liquid temperature = " << temperature_liquid << ",  ";
	cout << "solid temperature = " << temperature_solid << ",  ";
	cout << "liquid lamda = " << lamda_liquid << ",  ";
	cout << "density = " << pho << endl;
}