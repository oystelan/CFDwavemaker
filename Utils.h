#ifndef Utils_H
#define Utils_H

#include "irregular.h"

class Grid {
public:
	// Surface grid
	double* UX;
	double* UY;
	double* UZ;
	// lower grid used in the depth, where less exciting things are happening
	double* UXL;
	double* UYL;
	double* UZL;
	// Surface grid (2D)
	double* ETA;

	int NX, NY, NZ;
	int NXL, NYL, NZL;

	double domainsize[7];
	int initialized;
	int initsurf = 0;
	int initkin = 0;
	double dx, dy, dz, dxl, dyl, dzl;


	void initialize_kinematics(Irregular* irregular, double tpt);
	void initialize_surface_elevation(Irregular* irregular, double tpt);

	// Grid interpolation functions
	double trilinear_interpolation(double* VAR, double xpt, double ypt, double zpt);
	double trilinear_interpolationL(double* VAR, double xpt, double ypt, double zpt);
	double bilinear_interpolation(double* VAR, double xpt, double ypt);

};




#endif
