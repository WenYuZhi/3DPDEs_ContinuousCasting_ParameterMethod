#include <Mesh.h>
Mesh::Mesh(int nx, int ny, int nz, int tnpts, float lx, float ly, float lz, float tf)
{
	this->nx = nx;
	this->ny = ny;
	this->nz = nz;
	this->tnpts = tnpts;
	this->lx = lx;
	this->ly = ly;
	this->lz = lz;
	this->tf = tf;
	this->dx = this->lx / (this->nx - 1);
	this->dy = this->ly / (this->ny - 1);
	this->dz = this->lz / (this->nz - 1);
	this->tao = this->tf / (this->tnpts - 1);
}