#ifndef MESH_H
#define MESH_H

class Mesh
{
    public:
		int nx, ny, nz, tnpts;
		float dx, dy, dz, tao;
		float lx, ly, lz, tf;
		Mesh(int, int, int, int, float, float, float, float);
};
#endif
