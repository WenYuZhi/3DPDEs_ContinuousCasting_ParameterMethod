#include <Temperature.h>
#include <iostream>
#include <assert.h>
using namespace std;
Temperature::Temperature(const Mesh & continuous_casting_mesh, float vcast, const Continuous_Caster & casterone, Steel & steel)
{
	mCasterOne = &casterone;
	mesh = &continuous_casting_mesh;
	this->steel = &steel;
	T_New = new float[mesh->nx * mesh->ny * mesh->nz]();
	T_Last = new float[mesh->nx * mesh->ny * mesh->nz]();
	T_Surface = new float[mesh->ny]();
	meantemperature = new float[mCasterOne->coolsection]();
	this->vcast = vcast;
	tstep = 0;
	disout = true;
}

Temperature::Temperature(const Temperature & m_SteelTemperature)
{
	mCasterOne = m_SteelTemperature.mCasterOne;
	steel = m_SteelTemperature.steel;
	mesh = m_SteelTemperature.mesh;
	
	T_New = new float[mesh->nx * mesh->ny * mesh->nz]();
	T_Last = new float[mesh->nx * mesh->ny * mesh->nz]();
	for (int j = 0; j < mesh->ny; j++)
		for (int i = 0; i < mesh->nx; i++)
			for (int k = 0; k < mesh->nz; k++)
			{
				T_Last[mesh->nx * mesh->nz * j + mesh->nz * i + k] = m_SteelTemperature.T_Last[mesh->nx * mesh->nz * j + mesh->nz * i + k];
				T_New[mesh->nx * mesh->nz * j + mesh->nz * i + k] = m_SteelTemperature.T_New[mesh->nx * mesh->nz * j + mesh->nz * i + k];
			}
	T_Surface = new float[mesh->ny]();
	for (int j = 0; j < mesh->ny; j++)
		T_Surface[j] = m_SteelTemperature.T_Surface[j];

	meantemperature = new float[mCasterOne->coolsection]();
	for (int i = 0; i < mCasterOne->coolsection; i++)
		meantemperature[i] = m_SteelTemperature.meantemperature[i];

	vcast = m_SteelTemperature.vcast;
	tstep = m_SteelTemperature.tstep;
	disout = m_SteelTemperature.disout;
	h = m_SteelTemperature.h;
}

Temperature::~Temperature()
{
	delete[] T_New;
	delete[] T_Last;
	delete[] T_Surface;
	delete[] meantemperature;
}

void Temperature::setvcast(const int n_regression, const float* vcast_n_regression, vector<int> vcast_variation_tnpts)
{
	for (int i = 0; i < n_regression; i++)
		if (tstep >= vcast_variation_tnpts[i] && tstep < vcast_variation_tnpts[i + 1])
	        this->vcast = vcast_n_regression[i];
}

void Temperature::setcastingtemperature(float castingtemperature)
{
	this->castingtemperature = castingtemperature;
}

void Temperature::differencecalculation3d(const float *h_init, const int m_predictstep)
{
	float a, Tw = 30.0, T_Up, T_Down, T_Right, T_Left, T_Forw, T_Back, T_Middle;
	int nx = mesh->nx;
	int ny = mesh->ny;
	int nz = mesh->nz;
	float tao = mesh->tao;
	float dx = mesh->dx;
	float dy = mesh->dy;
	float dz = mesh->dz;

	for (int p = 0; p < m_predictstep; p++)
	{
		if (disout)
		{
			for (int j = 0; j < ny; j++)
			{
				h = this->boundarycondition_with_roll(h_init, j);
				for (int i = 0; i < nx; i++)
					for (int m = 0; m < nz; m++)
					{
						steel->physical_parameters(T_Last[nx * nz * j + nz * i + m]);
						a = steel->lamda / (steel->pho * steel->ce);
						if (j == 0 && i != 0 && i != (nx - 1) && m != 0 && m != (nz - 1)) //1
						{
							T_New[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == 0 && i == 0 && m != 0 && m != (nz - 1)) //2
						{
							T_New[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == 0 && i == (nx - 1) && m != 0 && m != (nz - 1))//3
						{
							T_New[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == 0 && i != 0 && i != (nx - 1) && m == 0) //4
						{
							T_New[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == 0 && i != 0 && i != (nx - 1) && m == (nz - 1)) //5
						{
							T_New[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == 0 && i == 0 && m == 0)  //6
						{
							T_New[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == 0 && i == 0 && m == (nz - 1))  //7
						{
							T_New[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == 0 && i == (nx - 1) && m == 0)  //8
						{
							T_New[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == 0 && i == (nx - 1) && m == (nz - 1)) //9
						{
							T_New[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m != 0 && m != (nz - 1)) //10
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
							T_Back = T_Last[nx * nz * j + nz * i + m - 1];
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == 0 && m != 0 && m != (nz - 1)) //11
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_Last[nx * nz * j + nz * (i + 1) + m];
							T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
							T_Back = T_Last[nx * nz * j + nz * i + m - 1];
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == (nx - 1) && m != 0 && m != (nz - 1)) //12
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i - 1) + m];
							T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
							T_Back = T_Last[nx * nz * j + nz * i + m - 1];
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m == 0)  //13
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
							T_Back = T_Last[nx * nz * j + nz * i + m + 1];
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m == (nz - 1))  //14
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m - 1];
							T_Back = T_Last[nx * nz * j + nz * i + m - 1];
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == 0 && m == 0)  //15
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_Last[nx * nz * j + nz * (i + 1) + m];
							T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
							T_Back = T_Last[nx * nz * j + nz * i + m + 1];
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == 0 && m == (nz - 1))  //16
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_Last[nx * nz * j + nz * (i + 1) + m];
							T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m - 1];
							T_Back = T_Last[nx * nz * j + nz * i + m - 1];
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == (nx - 1) && m == 0)  //17
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i - 1) + m];
							T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
							T_Back = T_Last[nx * nz * j + nz * i + m + 1];
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == (nx - 1) && m == (nz - 1))  //18
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i - 1) + m];
							T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m - 1];
							T_Back = T_Last[nx * nz * j + nz * i + m - 1];
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i != 0 && i != (nx - 1) && m == 0)  //19
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
							T_Back = T_Last[nx * nz * j + nz * i + m + 1] - 2 * dz * h * (T_Middle - Tw) / steel->lamda;
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i != 0 && i != (nx - 1) && m == (nz - 1))  //20
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m - 1] - 2 * dz * h * (T_Middle - Tw) / steel->lamda;
							T_Back = T_Last[nx * nz * j + nz * i + m - 1];
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == 0 && m == 0) //21
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_Last[nx * nz * j + nz * (i + 1) + m];
							T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
							T_Back = T_Last[nx * nz * j + nz * i + m + 1];
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m == 0)  //22
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i - 1) + m];
							T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
							T_Back = T_Last[nx * nz * j + nz * i + m + 1];
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == 0 && m == (nz - 1)) //23
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_Last[nx * nz * j + nz * (i + 1) + m];
							T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m - 1];
							T_Back = T_Last[nx * nz * j + nz * i + m - 1];
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m == (nz - 1)) //24
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i - 1) + m];
							T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m - 1];
							T_Back = T_Last[nx * nz * j + nz * i + m - 1];
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == 0 && m != 0 && m != (nz - 1))  //25
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_Last[nx * nz * j + nz * (i + 1) + m] - 2 * dx * h * (T_Middle - Tw) / steel->lamda;
							T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
							T_Back = T_Last[nx * nz * j + nz * i + m - 1];
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m != 0 && m != (nz - 1)) //26
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i - 1) + m] - 2 * dx * h * (T_Middle - Tw) / steel->lamda;
							T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
							T_Back = T_Last[nx * nz * j + nz * i + m - 1];
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else  //27
						{
							T_Middle = T_Last[nx * nz * j + nz * i + m];
							T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
							T_Back = T_Last[nx * nz * j + nz * i + m - 1];
							T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}
					}
			}
			for (int k = 0; k < ny; k++)
				T_Surface[k] = T_New[nx * nz * k + nz * int((nx - 1) / 2) + nz - 1];
		}

		else
		{
			for (int j = 0; j < ny; j++)
			{
				h = this->boundarycondition_with_roll(h_init, j);
				for (int i = 0; i < nx; i++)
					for (int m = 0; m < nz; m++)
					{
						steel->physical_parameters(T_New[nx * nz * j + nz * i + m]);
						a = steel->lamda / (steel->pho * steel->ce);
						if (j == 0 && i != 0 && i != (nx - 1) && m != 0 && m != (nz - 1)) //1
						{
							T_Last[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == 0 && i == 0 && m != 0 && m != (nz - 1)) //2
						{
							T_Last[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == 0 && i == (nx - 1) && m != 0 && m != (nz - 1))//3
						{
							T_Last[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == 0 && i != 0 && i != (nx - 1) && m == 0) //4
						{
							T_Last[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == 0 && i != 0 && i != (nx - 1) && m == (nz - 1)) //5
						{
							T_Last[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == 0 && i == 0 && m == 0)  //6
						{
							T_Last[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == 0 && i == 0 && m == (nz - 1))  //7
						{
							T_Last[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == 0 && i == (nx - 1) && m == 0)  //8
						{
							T_Last[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == 0 && i == (nx - 1) && m == (nz - 1)) //9
						{
							T_Last[nx * nz * j + nz * i + m] = castingtemperature;
						}

						else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m != 0 && m != (nz - 1)) //10
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == 0 && m != 0 && m != (nz - 1)) //11
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == (nx - 1) && m != 0 && m != (nz - 1)) //12
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m == 0)  //13
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m + 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m == (nz - 1))  //14
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m - 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == 0 && m == 0)  //15
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m + 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == 0 && m == (nz - 1))  //16
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m - 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == (nx - 1) && m == 0)  //17
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m + 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == (nx - 1) && m == (nz - 1))  //18
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m - 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i != 0 && i != (nx - 1) && m == 0)  //19
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m + 1] - 2 * dz * h * (T_Middle - Tw) / steel->lamda;
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i != 0 && i != (nx - 1) && m == (nz - 1))  //20
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m - 1] - 2 * dz * h * (T_Middle - Tw) / steel->lamda;
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == 0 && m == 0) //21
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m + 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m == 0)  //22
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m + 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == 0 && m == (nz - 1)) //23
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m - 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m == (nz - 1)) //24
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m - 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == 0 && m != 0 && m != (nz - 1))  //25
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i + 1) + m] - 2 * dx * h * (T_Middle - Tw) / steel->lamda;
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m != 0 && m != (nz - 1)) //26
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i - 1) + m] - 2 * dx * h * (T_Middle - Tw) / steel->lamda;
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else  //27
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}
					}
			}
			for (int k = 0; k < ny; k++)
				T_Surface[k] = T_Last[nx * nz * k + nz * int((nx - 1) / 2) + nz - 1];
		}
		disout = !disout;
		tstep++;
	}
}

float Temperature::boundarycondition3d(const float *h_init, int j)
{
	float yposition = j * mesh->dy, h = 0.0f;
	for (int i = 0; i < mCasterOne->section; i++)
		if (yposition >= *(mCasterOne->ccml + i) && yposition <= *(mCasterOne->ccml + i + 1))
			h = *(h_init + i);
	return h;
}

float Temperature::boundarycondition_with_roll(const float *h_init, int j)
{
	float yposition = j * mesh->dy, h_rad = 100.0f, h_roll = 0.0f, h_spray = 0.0f;
	for (int i = 0; i < mCasterOne->section; i++)
	{
		if (yposition >= *(mCasterOne->ccml + i) && yposition <= *(mCasterOne->ccml + i + 1))
			h_spray = *(h_init + i);
	}

	for (int i = 0; i < mCasterOne->n_roll; i++)
		if (fabs(yposition - mCasterOne->roll_position[i]) < mCasterOne->length_roll_contact)
		{
			h_roll = 0.05f / ((1.0f - 0.05f) * mCasterOne->length_roll_contact) * (h_rad + h_spray) * (mCasterOne->roll_position[i + 1] - mCasterOne->roll_position[i] - mCasterOne->length_roll_contact);
			return h_roll;
		}

	return h_spray;
}

void Temperature::initcondition3d(float castingtemperature)
{
	tstep = 0;
	for (int j = 0; j < mesh->ny; j++)
		for (int i = 0; i < mesh->nx; i++)
			for (int k = 0; k < mesh->nz; k++)
			{
				T_Last[mesh->nx * mesh->nz * j + mesh->nz * i + k] = castingtemperature;
				T_New[mesh->nx * mesh->nz * j + mesh->nz * i + k] = castingtemperature;
			}
	disout = true;
}

void Temperature::initcondition3d(float *castingtemperature)
{
	tstep = 0;
	for (int j = 0; j < mesh->ny; j++)
		for (int i = 0; i < mesh->nx; i++)
			for (int k = 0; k < mesh->nz; k++)
			{
				T_Last[mesh->nx * mesh->nz * j + mesh->nz * i + k] = castingtemperature[mesh->nx * mesh->nz * j + mesh->nz * i + k];
				T_New[mesh->nx * mesh->nz * j + mesh->nz * i + k] = castingtemperature[mesh->nx * mesh->nz * j + mesh->nz * i + k];
			}
	disout = true;
}


void Temperature::computemeantemperature3d()
{
	float y;
	int count = 0;
	for (int i = 0; i < mCasterOne->coolsection; i++)
	{
		meantemperature[i] = 0.0;
		for (int j = 0; j < mesh->ny; j++)
		{
			y = j * mesh->dy;
			if (y > *((mCasterOne->ccml) + i + mCasterOne->moldsection) && y <= *((mCasterOne->ccml) + i + 1 + +mCasterOne->moldsection))
			{
				meantemperature[i] += T_Surface[j];
				count++;
			}
		}
		assert(count != 0);
		meantemperature[i] = meantemperature[i] / count;
		count = 0;
	}

}

void Temperature::print3d(int measurednumb)
{
	if (tstep == 0)
	{
		cout << "lx = " << mesh->lx << ", " << "nx = " << mesh->nx << ", ";
		cout << "ly = " << mesh->ly << ", " << "ny = " << mesh->ny << ", ";
		cout << "lz = " << mesh->lz << ", " << "nz = " << mesh->nz << ", ";
		cout << "casting speed = " << vcast << ", " << endl;
		cout << "dx = " << mesh->dx << ", ";
		cout << "dy = " << mesh->dy << ", ";
		cout << "dz = " << mesh->dz << ", ";
		cout << "time step = " << mesh->tao << ", " << endl;
	}
}

void Temperature::print3d()
{
	if (tstep == 0)
	{
		cout << "lx = " << mesh->lx << ", " << "nx = " << mesh->nx << ", ";
		cout << "ly = " << mesh->ly << ", " << "ny = " << mesh->ny << ", ";
		cout << "lz = " << mesh->lz << ", " << "nz = " << mesh->nz << ", ";
		cout << "casting speed = " << vcast << ", " << endl;
		cout << "dx = " << mesh->dx << ", ";
		cout << "dy = " << mesh->dy << ", ";
		cout << "dz = " << mesh->dz << ", ";
		cout << "time step = " << mesh->tao << ", " << endl;
	}
	else
	{
		cout << "tstep = " << tstep << endl;
		cout << "meantemperature = " << endl;
		for (int i = 0; i < mCasterOne->coolsection; i++)
			cout << meantemperature[i] << ", ";
		cout << endl;
	}
}

void Temperature::operator=(const Temperature & m_SteelTemperature)
{
	mCasterOne = m_SteelTemperature.mCasterOne;
	steel = m_SteelTemperature.steel;
	mesh = m_SteelTemperature.mesh;
	vcast = m_SteelTemperature.vcast;
	castingtemperature = m_SteelTemperature.castingtemperature;
	delete[] T_New;
	delete[] T_Last;
	delete[] T_Surface;
	delete[] meantemperature;

	T_New = new float[mesh->nx *mesh->ny *mesh->nz];
	T_Last = new float[mesh->nx *mesh->ny *mesh->nz];
	for (int j = 0; j < mesh->ny; j++)
		for (int i = 0; i < mesh->nx; i++)
			for (int k = 0; k < mesh->nz; k++)
			{
				T_Last[mesh->nx * mesh->nz * j + mesh->nz * i + k] = m_SteelTemperature.T_Last[mesh->nx * mesh->nz * j + mesh->nz * i + k];
				T_New[mesh->nx * mesh->nz * j + mesh->nz * i + k] = m_SteelTemperature.T_New[mesh->nx * mesh->nz * j + mesh->nz * i + k];
			}
	T_Surface = new float[mesh->ny];
	for (int j = 0; j < mesh->ny; j++)
		T_Surface[j] = m_SteelTemperature.T_Surface[j];

	meantemperature = new float[mCasterOne->coolsection];
	for (int i = 0; i < mCasterOne->coolsection; i++)
		meantemperature[i] = m_SteelTemperature.meantemperature[i];

	vcast = m_SteelTemperature.vcast;
	tstep = m_SteelTemperature.tstep;
	disout = m_SteelTemperature.disout;
}