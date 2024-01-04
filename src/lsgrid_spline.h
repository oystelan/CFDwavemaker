#ifndef PI
#define PI 3.1415926535897
#endif

#ifndef SGrid_H
#define SGrid_H
#include <vector>
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
		delete[] UX, UY, UZ, ETA, IGNORE;
	};

	// Primary field variables
	double* UX, * UY, * UZ; // 4D dimensions, (time, x, y, z)	
	double* ETA;// Surface grid (3D) (time, x, y)
	// Secondary field variables (derivatives)
	double* ETAdx,* ETAdy,* ETAdt;
	double* Udt, * Udx,* Udy,* Uds;
	double* Vdt, * Vdx,* Vdy,* Vds;
	double* Wdt, * Wdx,* Wdy,* Wds;

	int* IGNORE; // matrix with cells to ignore when updating surface elevation
	int nx, ny, nl;
	double domain[4] = {};
	double domain_ignore[4] = {largeval, -largeval, largeval, -largeval};
	bool ignore_at_init = false;
	int tstep = 0;

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

	double sigmaS2cart(double s, double wave_elev, double depth);
	double cart2sigmaS(double zpt, double wave_elev, double depth);

	void square_vals(double* C, double* data, int nxp, int nyp, int tid);
	void cube_vals(double* C, double* data, int nxp, int nyp, int nlp, int tid);
	double spline_interp_velo(double* U, double* Udt, double* Udx, double* Udy, double* Uds, int nxp, int nyp, int nsp0, int nsp1, double xd, double yd, double sd0, double sd1, double td, int tid1, int tid2);

	
	void write_vtk(bool endtime);
	void allocate();

	// Second order theory
	void initialize_kinematics(Irregular& irregular);
	void initialize_kinematics_with_ignore(Irregular& irregular);
	void update(Irregular& irregular, double t_target);

/*
	// SWD
#if defined(SWD_enable)
	void initialize_surface_elevation_with_ignore(SpectralWaveData* swd, double t_target);
	void initialize_surface_elevation(SpectralWaveData* swd, double t_target);
	void initialize_kinematics_with_ignore(SpectralWaveData* swd);
	void initialize_kinematics(SpectralWaveData* swd);
	void update(SpectralWaveData* swd, double t_target);
#endif
*/

	void set_ignore();
	// gradient update functions
	void update_gradient_dt(double* data, double* graddt);
	void update_gradient_eta_dt(double* data, double* graddt);
	void update_gradient_eta_dxdy(double* eta, double* gradx,double* grady);
	void update_gradient_dxdydz(double* data, double* gradx, double* grady, double* gradz);

	// Grid interpolation functions
	std::vector<double> get_kinematics_at_point(double tpt, double xpt, double ypt, double zpt, double h);


	// routines for extracting kinematics at walls
	bool CheckTime(double tpt);
	bool CheckBounds();

	void update_bounds(double xpt, double ypt);
	
	// export_functions to validate kinematics
	void export_vtu(FILE* fp, bool last);

};

#endif