#ifndef Utils_H
#define Utils_H

#include "irregular.h"

class Grid {
public:
	Grid() {

	};
	~Grid() {
		delete[] UX, UY, UZ, UXL, UYL, UZL, ETA;
		delete[] UX0, VX0, WX0, UX1, VX1, WX1;
		delete[] UY0, VY0, WY0, UY1, VY1, WY1;
		delete[] ETAX0, ETAX1, ETAY0, ETAY1;
	};
	// Volume grid, used for fast initialization
	double * UX, * UY, * UZ, * UXL, * UYL, * UZL; // 3D grids	
	double* ETA;// Surface grid (2D)

	
	// Wall grids, used for faster running of calculation demanding kinematics codes
	double* UX0, * VX0, * WX0, * UX1, * VX1, * WX1; // X wall grid (2D)
	double* UY0, * VY0, * WY0, * UY1, * VY1, * WY1; // Y wall grid (2D)
	double* ETAX0, * ETAX1, * ETAY0, * ETAY1; // surface grids at wall (1D)

	int NX, NY, NZ;
	int NXL, NYL, NZL;

	double domainsize[7];
	double wallxsize[6];
	double wallysize[6];
	int wallx_nx, wallx_nz;
	int wally_ny, wally_nz;

	int initialized;
	int initsurf = 0;
	int initkin = 0;
	double dx, dy, dz, dxl, dyl, dzl;
	int numgrids;

	void initialize_kinematics(Irregular* irregular, double tpt);
	void initialize_surface_elevation(Irregular* irregular, double tpt);

	// Grid interpolation functions
	double trilinear_interpolation(double* VAR, double xpt, double ypt, double zpt);
	double trilinear_interpolationL(double* VAR, double xpt, double ypt, double zpt);
	double bilinear_interpolation(double* VAR, double xpt, double ypt);

};
// Ramp class contains various ramp functions
class Ramp {
private:
	double ramp1d(double x, double xstart, double xend, bool inv);
public:
	int timeramp, xramp, yramp; // logical operators for switching on and of rampfunctions
	double time_rampup_start, time_rampup_end, time_rampdown_start, time_rampdown_end;
	double x_rampup_start, x_rampup_end, x_rampdown_start, x_rampdown_end;
	double y_rampup_start, y_rampup_end, y_rampdown_start, y_rampdown_end;
	bool ramp_init = false;
	bool ramp_init_time_up = false;
	bool ramp_init_time_down = false;
	bool ramp_init_x_up = false;
	bool ramp_init_x_down = false;
	bool ramp_init_y_up = false;
	bool ramp_init_y_down = false;

	double ramp(double time, double x, double y);
};




#endif
