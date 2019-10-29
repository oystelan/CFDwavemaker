#ifndef Utils_H
#define Utils_H

#include "irregular.h"

class Grid {
public:
	// Volume grid, used for fast initialization
	double* UX, UY, UZ, UXL, UYL, UZL;

	// Wall grids, used for faster running of calculation demanding kinematics codes
	double* UX0, UY0, UZ0, UX1, UY1, UZ1;

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
