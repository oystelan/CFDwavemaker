#ifndef PI
#define PI 3.1415926535897
#endif

#ifndef SGrid_H
#define SGrid_H

#include "Irregular.h"


class sGrid {
public:
	sGrid() {
		// Set default values
		tan_a = 7. * PI / 18.;
		tan_b = 1.5;
		NL = 10;
	};
	~sGrid() {
		delete[] UX0, UY0, UZ0, ETA0, UX1, UY1, UZ1, IGNORE;
	};
	// Volume grid, used for fast initialization
	double* UX0, * UY0, * UZ0, * UX1, * UY1, * UZ1; // 3D grids	
	double* ETA0, * ETA1;// Surface grid (2D)
	int* IGNORE; // matrix with cells to ignore when updating surface elevation
	int NX, NY, NL;
	double domain[4] = {};
	double domain_ignore[4] = {};

	double water_depth;
	double dx, dy, ds;
	double t0, dt;
	int numgrids;
	double tan_a;
	double tan_b;

	double slayer(int layerno);
	double clayer(int layerno);

	// sGrid transform functions
	double z2s(double z, double wave_elev, double depth);
	double s2z(double s, double wave_elev, double depth);
	double s2tan(double s);
	double tan2s(double t);

	void initialize_kinematics(Irregular& irregular);
	void initialize_surface_elevation(Irregular& irregular);

	// Update kinematics data matrices
	void update(Irregular& irregular);

	void set_ignore(double* bounds);

	// Grid interpolation functions
	double trilinear_interpolation(double* VAR, double xpt, double ypt, double zpt);
	double bilinear_interpolation(double* VAR, double xpt, double ypt);

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