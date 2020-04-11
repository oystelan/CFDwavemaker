#ifndef SGrid_H
#define SGrid_H

#include "Irregular.h"

class sGrid {
public:
	sGrid() {

	};
	~sGrid() {
		delete[] UX0, UY0, UZ0, ETA0;
	};
	// Volume grid, used for fast initialization
	double* UX0, * UY0, * UZ0, * UX1, * UY1, * UZ1; // 3D grids	
	double* ETA0, * ETA0;// Surface grid (2D)
	int NX, NY, NZ;
	double domain_start[3], domain_end[3];

	double water_depth;
	bool initialized;
	int initsurf = 0;
	int initkin = 0;
	int nl = 10; // number of layers in z direction
	double dx, dy, dz;
	double t0, t1, dt;
	int numgrids;

	double slayer(int layerno);

	double clayer(int layerno);

	double z2s(double z, double wave_elev, double depth);

	double s2z(double s, double wave_elev, double depth);

	double s2tan(double s);

	double tan2s(double t);

	void initialize_kinematics(Irregular& irregular, double tpt);
	void initialize_surface_elevation(Irregular& irregular, double tpt);

	// Grid interpolation functions
	double trilinear_interpolation(double* VAR, double xpt, double ypt, double zpt, int _nx, int _ny, int _nz, double _dx, double _dy, double _dz, double* domain);
	double bilinear_interpolation(double* VAR, double xpt, double ypt, int _nx, int _ny, double _dx, double _dy, double* domain);

	// routines for extracting kinematics from initial grid
	double u(double xpt, double ypt, double zpt);
	double v(double xpt, double ypt, double zpt);
	double w(double xpt, double ypt, double zpt);
	double eta(double xpt, double ypt);

	// routines for extracting kinematics at walls
	bool CheckTime(double tpt);
	bool CheckBounds(double* bounds, double x, double y, double z);
};

#endif