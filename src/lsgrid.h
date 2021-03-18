#ifndef PI
#define PI 3.1415926535897
#endif

#ifndef SGrid_H
#define SGrid_H

#include "Irregular.h"
#if SWD_enable
#include "SpectralWaveData.h"
#endif

#define largeval 1.E12

class lsGrid {
public:
	lsGrid() {
		// Set default values
		tan_a = 7. * PI / 18.;
		tan_b = 1.5;
		nl = 10;
		dt = 0.5;
		swl = 0.;
	};
	~lsGrid() {
		if (dump_vtk) {
			write_vtk(true);
		}
		delete[] UX0, UY0, UZ0, ETA0, UX1, UY1, UZ1, IGNORE;
	};

	// Volume grid, used for fast initialization
	double* UX0, * UY0, * UZ0, * UX1, * UY1, * UZ1; // 3D grids	
	double* ETA0, * ETA1;// Surface grid (2D)
	int* IGNORE; // matrix with cells to ignore when updating surface elevation
	int nx, ny, nl;
	double domain[4] = {};
	double domain_ignore[4] = {largeval, -largeval, largeval, -largeval};
	bool ignore_at_init = false;
	int update_count = 0;

	bool disable_checkbounds = true;
	bool ignore_domain = false;
	bool init_only = false;

	std::string vtk_directory_path = "./";
	std::string vtk_prefix = "kinematics";
	std::string vtk_timelabel = "TimeValue";
	bool dump_vtk = false;

	double water_depth;
	double swl;
	double dx, dy, ds;
	double t0, dt;
	double tan_a;
	double tan_b;

	double bxmin, bxmax, bymin, bymax;

	double slayer(int layerno);
	double clayer(int layerno);

	// sGrid transform functions
	double z2s(double z, double wave_elev, double depth);
	double s2z(double s, double wave_elev, double depth);
	double s2tan(double s);
	double tan2s(double t);

	
	
	void write_vtk(bool endtime);
	void allocate();

	// Second order theory
	void initialize_kinematics(Irregular& irregular);
	void initialize_kinematics_with_ignore(Irregular& irregular);
	void initialize_surface_elevation(Irregular& irregular, double t_target);
	void initialize_surface_elevation_with_ignore(Irregular& irregular, double t_target);
	void update(Irregular& irregular, double t_target);

	// SWD
#if defined(SWD_enable)
	void initialize_surface_elevation_with_ignore(SpectralWaveData* swd, double t_target);
	void initialize_surface_elevation(SpectralWaveData* swd, double t_target);
	void initialize_kinematics_with_ignore(SpectralWaveData* swd);
	void initialize_kinematics(SpectralWaveData* swd);
	void update(SpectralWaveData* swd, double t_target);
#endif

	void set_ignore();

	// Grid interpolation functions
	double trilinear_interpolation(double* VAR0, double* VAR1, double tpt, double xpt, double ypt, double zpt);
	double trilinear_interpolation2(double* VAR0, double* VAR1, double tpt, double xpt, double ypt, double zpt);
	double bilinear_interpolation(double* VAR0, double* VAR1, double tpt, double xpt, double ypt);

	double bilinear_interpolation2(double* VAR0, double* VAR1, double tpt, double xpt, double ypt);

	// routines for extracting kinematics from initial grid
	double u(double tpt, double xpt, double ypt, double zpt);
	double v(double tpt, double xpt, double ypt, double zpt);
	double w(double tpt, double xpt, double ypt, double zpt);
	double eta(double tpt, double xpt, double ypt);

	// routines for extracting kinematics at walls
	bool CheckTime(double tpt);
	bool CheckBounds();

	void update_bounds(double xpt, double ypt);
	
	// export_functions to validate kinematics
	void export_vtu(FILE* fp, bool last);

};

#endif