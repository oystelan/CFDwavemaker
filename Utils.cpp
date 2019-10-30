
#include "Utils.h"
#include "omp.h"
#include <algorithm>
#include <iostream>


// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 3
void Grid::initialize_kinematics(Irregular *irregular, double tpt) {
	t0 = tpt;
	t1 = tpt;
	// Allocating memory for storage of surface elevation and velocities
	UX = new double[NX * NY * NZ];
	UY = new double[NX * NY * NZ];
	UZ = new double[NX * NY * NZ];

	UXL = new double[NXL * NYL * NZL];
	UYL = new double[NXL * NYL * NZL];
	UZL = new double[NXL * NYL * NZL];

	std::cout << "Memory allocation successful for storage of kinematics." << std::endl;

	dx = (domain_end[0] - domain_start[0]) / double(NX - 1);
	dy = (domain_end[1] - domain_start[1]) / double(NY - 1);
	dz = (domain_end[2] - domain_start[2]) / double(NZ - 1);

	double dd = omp_get_wtime();

	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());

#pragma omp parallel // start parallel initialization
	{
#pragma omp master
		std::cout << "Number of available threads: " << omp_get_num_threads() << std::endl;

		double xpt, ypt, zpt;
		double eta_temp;

		// Main grid
#pragma omp for
		for (int i = 0; i < NX; i++) {
			xpt = domain_start[0] + dx * i;
			for (int j = 0; j < NY; j++) {
				ypt = domain_start[1] + dy * j;
				eta_temp = irregular->eta(tpt, xpt, ypt);

				double Ux0 = irregular->u1(tpt, xpt, ypt, 0.0) + irregular->u2(tpt, xpt, ypt, 0.0);
				double Uy0 = irregular->v1(tpt, xpt, ypt, 0.0) + irregular->v2(tpt, xpt, ypt, 0.0);
				double Uz0 = irregular->w1(tpt, xpt, ypt, 0.0) + irregular->w2(tpt, xpt, ypt, 0.0);

				double PHI_dxdz = irregular->phi1_dxdz(tpt, xpt, ypt);
				double PHI_dydz = irregular->phi1_dydz(tpt, xpt, ypt);
				double PHI_dzdz = irregular->phi1_dzdz(tpt, xpt, ypt);

				for (int m = 0; m < NZ; m++) {
					zpt = domain_start[2] + dz * m;
					if (zpt > (eta_temp + dz)) {
						UX[i * NY * NZ + j * NZ + m] = 0.0;
						UY[i * NY * NZ + j * NZ + m] = 0.0;
						UZ[i * NY * NZ + j * NZ + m] = 0.0;
					}
					else if (zpt > 0.) {
						UX[i * NY * NZ + j * NZ + m] = Ux0 + PHI_dxdz * zpt;
						UY[i * NY * NZ + j * NZ + m] = Uy0 + PHI_dydz * zpt;
						UZ[i * NY * NZ + j * NZ + m] = Uz0 + PHI_dzdz * zpt;
					}
					else {
						UX[i * NY * NZ + j * NZ + m] = irregular->u1(tpt, xpt, ypt, zpt) + irregular->u2(tpt, xpt, ypt, zpt);
						UY[i * NY * NZ + j * NZ + m] = irregular->v1(tpt, xpt, ypt, zpt) + irregular->v2(tpt, xpt, ypt, zpt);
						UZ[i * NY * NZ + j * NZ + m] = irregular->w1(tpt, xpt, ypt, zpt) + irregular->w2(tpt, xpt, ypt, zpt);
					}
					/*UX[i*NY*NZ + j*NZ + m] = uu(tpt, xpt, ypt, zpt);
					UY[i*NY*NZ + j*NZ + m] = vv(tpt, xpt, ypt, zpt);
					UZ[i*NY*NZ + j*NZ + m] = ww(tpt, xpt, ypt, zpt);*/
				}
			}
		}
	} // End parallel initialization

	std::cout << "Generation of upper domain kinematics data completed. ";
	dd = omp_get_wtime() - dd;
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;

	dd = omp_get_wtime();
#pragma omp parallel // start parallel initialization
	{
		double xpt, ypt, zpt;
		// Secondary grid (coarse res at depth)
		dxl = (domain_end_L[0] - domain_start_L[0]) / double(NXL - 1);
		dyl = (domain_end_L[1] - domain_start_L[1]) / double(NYL - 1);
		dzl = (domain_end_L[2] - domain_start_L[2]) / double(NZL - 1);
#pragma omp for
		for (int i = 0; i < NXL; i++) {
			xpt = domain_start_L[0] + dxl * i;
			for (int j = 0; j < NYL; j++) {
				ypt = domain_start_L[1] + dyl * j;
				for (int m = 0; m < NZL; m++) {
					zpt = domain_start_L[2] + dzl * m;
					UXL[i * NYL * NZL + j * NZL + m] = irregular->u1(tpt, xpt, ypt, zpt) + irregular->u2(tpt, xpt, ypt, zpt);
					UYL[i * NYL * NZL + j * NZL + m] = irregular->v1(tpt, xpt, ypt, zpt) + irregular->v2(tpt, xpt, ypt, zpt);
					UZL[i * NYL * NZL + j * NZL + m] = irregular->w1(tpt, xpt, ypt, zpt) + irregular->w2(tpt, xpt, ypt, zpt);
				}
			}
		}
	} // End parallel
	std::cout << "Generation of lower domain kinematics data completed. ";
	dd = omp_get_wtime() - dd;
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
	std::cout << "Interpolation can commence..." << std::endl;
	initkin = 1;
}


// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 3
void Grid::initialize_surface_elevation(Irregular* irregular, double tpt) {

	// Allocating memory for storage of surface elevation and velocities
	ETA = new double[NX * NY];

	std::cout << "Memory allocation successful for Surface elevation storage." << std::endl;

	dx = (domain_end[0] - domain_start[0]) / double(NX - 1);
	dy = (domain_end[1] - domain_start[1]) / double(NY - 1);

	double dd = omp_get_wtime();
	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());

#pragma omp parallel
	{
		double xpt, ypt;
		// Main grid
#pragma omp for
		for (int i = 0; i < NX; i++) {
			xpt = domain_start[0] + dx * i;
			for (int j = 0; j < NY; j++) {
				ypt = domain_start[1] + dy * j;
				ETA[i * NY + j] = irregular->eta1(tpt, xpt, ypt) + irregular->eta2(tpt, xpt, ypt);
			}
		}
	}
	dd = omp_get_wtime() - dd;

	std::cout << "Surface Elevation generated successfully. ";
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
	initsurf = 1;
}

/* Function for trilinear interpolation on a cartesian evenly spaced mesh*/
double Grid::trilinear_interpolation(double* VAR, double xpt, double ypt, double zpt, int _nx, int _ny, int _nz, double _dx, double _dy, double _dz, double *domain) {
	double nxp = std::min(double(_nx), std::max(0., (xpt - domain[0]) / _dx));
	double nyp = std::min(double(_ny), std::max(0., (ypt - domain[1]) / _dy));
	double nzp = std::min(double(_nz), std::max(0., (zpt - domain[2]) / _dz));

	double C000 = VAR[int(floor(nxp) * _ny * _nz + floor(nyp) * _nz + floor(nzp))];
	double C001 = VAR[int(floor(nxp) * _ny * _nz + floor(nyp) * _nz + ceil(nzp))];
	double C010 = VAR[int(floor(nxp) * _ny * _nz  + ceil(nyp) * _nz + floor(nzp))];
	double C011 = VAR[int(floor(nxp) * _ny * _nz + ceil(nyp) * _nz + ceil(nzp))];
	double C100 = VAR[int(ceil(nxp) * _ny * _nz + floor(nyp) * _nz + floor(nzp))];
	double C101 = VAR[int(ceil(nxp) * _ny * _nz + floor(nyp) * _nz + ceil(nzp))];
	double C110 = VAR[int(ceil(nxp) * _ny * _nz + ceil(nyp) * _nz + floor(nzp))];
	double C111 = VAR[int(ceil(nxp) * _ny * _nz + ceil(nyp) * _nz + ceil(nzp))];
	double xd = nxp - floor(nxp);
	double yd = nyp - floor(nyp);
	double zd = nzp - floor(nzp);

	double C00 = C000 * (1. - xd) + C100 * xd;
	double C01 = C001 * (1. - xd) + C101 * xd;
	double C10 = C010 * (1. - xd) + C110 * xd;
	double C11 = C011 * (1. - xd) + C111 * xd;

	double C0 = C00 * (1. - yd) + C10 * yd;
	double C1 = C01 * (1. - yd) + C11 * yd;

	return C0 * (1. - zd) + C1 * zd;
}


/* Function for trilinear interpolation on a cartesian evenly spaced mesh on the lower part of the domain*/
/*double Grid::trilinear_interpolationL(double* VAR, double xpt, double ypt, double zpt) {

	double nxp = std::min(double(NXL), std::max(0., (xpt - domainsize[0]) / dxl));
	double nyp = std::min(double(NYL), std::max(0., (ypt - domainsize[2]) / dyl));
	double nzp = std::min(double(NZL), std::max(0., (zpt - domainsize[4]) / dzl));

	double C000 = VAR[int(floor(nxp) * NYL * NZL + floor(nyp) * NZL + floor(nzp))];
	double C001 = VAR[int(floor(nxp) * NYL * NZL + floor(nyp) * NZL + ceil(nzp))];
	double C010 = VAR[int(floor(nxp) * NYL * NZL + ceil(nyp) * NZL + floor(nzp))];
	double C011 = VAR[int(floor(nxp) * NYL * NZL + ceil(nyp) * NZL + ceil(nzp))];
	double C100 = VAR[int(ceil(nxp) * NYL * NZL + floor(nyp) * NZL + floor(nzp))];
	double C101 = VAR[int(ceil(nxp) * NYL * NZL + floor(nyp) * NZL + ceil(nzp))];
	double C110 = VAR[int(ceil(nxp) * NYL * NZL + ceil(nyp) * NZL + floor(nzp))];
	double C111 = VAR[int(ceil(nxp) * NYL * NZL + ceil(nyp) * NZL + ceil(nzp))];
	double xd = nxp - floor(nxp);
	double yd = nyp - floor(nyp);
	double zd = nzp - floor(nzp);

	double C00 = C000 * (1. - xd) + C100 * xd;
	double C01 = C001 * (1. - xd) + C101 * xd;
	double C10 = C010 * (1. - xd) + C110 * xd;
	double C11 = C011 * (1. - xd) + C111 * xd;

	double C0 = C00 * (1. - yd) + C10 * yd;
	double C1 = C01 * (1. - yd) + C11 * yd;

	return C0 * (1. - zd) + C1 * zd;
}*/

/* bilinear interpolation function used to interpolate surface values on a regular evenly spaced grid*/
double Grid::bilinear_interpolation(double* VAR, double xpt, double ypt, int _nx, int _ny, double _dx, double _dy, double *domain) {

	double nxp = std::min(double(_nx), std::max(0., (xpt - domain[0]) / _dx));
	double nyp = std::min(double(_ny), std::max(0., (ypt - domain[1]) / _dy));

	double C00 = VAR[int(floor(nxp) * _ny + floor(nyp))];
	double C01 = VAR[int(floor(nxp) * _ny + ceil(nyp))];
	double C10 = VAR[int(ceil(nxp) * _ny + floor(nyp))];
	double C11 = VAR[int(ceil(nxp) * _ny + ceil(nyp))];

	double xd = nxp - floor(nxp);
	double yd = nyp - floor(nyp);

	double C0 = C00 * (1. - xd) + C10 * xd;
	double C1 = C01 * (1. - xd) + C11 * xd;

	return C0 * (1. - yd) + C1 * yd;
}


double Grid::u(double xpt, double ypt, double zpt) {
	if (zpt < domain_start[2]) {
		return trilinear_interpolation(UXL, xpt, ypt, zpt, NXL, NYL, NZL, dxl, dyl, dzl, domain_start_L);
	}
	else {
		return trilinear_interpolation(UX, xpt, ypt, zpt, NX, NY, NZL, dx, dy, dz, domain_start);
	}
}

double Grid::v(double xpt, double ypt, double zpt) {
	if (zpt < domain_start[2]) {
		return trilinear_interpolation(UYL, xpt, ypt, zpt, NXL, NYL, NZL, dxl, dyl, dzl, domain_start_L);
	}
	else {
		return trilinear_interpolation(UY, xpt, ypt, zpt, NX, NY, NZL, dx, dy, dz, domain_start);
	}
}

double Grid::w(double xpt, double ypt, double zpt) {
	if (zpt < domain_start[2]) {
		return trilinear_interpolation(UZL, xpt, ypt, zpt, NXL, NYL, NZL, dxl, dyl, dzl, domain_start_L);
	}
	else {
		return trilinear_interpolation(UZ, xpt, ypt, zpt, NX, NY, NZL, dx, dy, dz, domain_start);
	}
}

double Grid::eta(double xpt, double ypt) {
	return bilinear_interpolation(UZ, xpt, ypt, NX, NY, dx, dy, domain_start);
}




// -------------------------------------------------------------------------------------------------
// ramp class function
// -------------------------------------------------------------------------------------------------

//Define some useful functions
/* Rampfunction */
// NB: Not yet implemented inverse ramp
double Ramp::ramp1d(double x, double xstart, double xend, bool inv) {

	if (inv){
		return std::max(std::min(1. - ((x - xstart) / (xend - xstart)),1.), 0.);
	}
	else {
		return std::max(std::min(((x - xstart) / (xend - xstart)), 1.), 0.);
	}
}

bool comp(double a, double b)
{
	return (a < b);
}

// Todo: implement ramps
double Ramp::ramp(double t, double x, double y) {
	if (ramp_init) {
		double ramps[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
		if (ramp_init_time_up) {
			ramps[0] = ramp1d(t, time_rampup_start, time_rampup_end, false);
		}
		if (ramp_init_time_down) {
			ramps[1] = ramp1d(t, time_rampdown_start, time_rampdown_end, true);
		}
		if (ramp_init_x_up) {
			ramps[2] = ramp1d(x, x_rampup_start, x_rampup_end, false);
		}
		if (ramp_init_x_down) {
			ramps[3] = ramp1d(x, x_rampdown_start, x_rampdown_end, true);
		}
		if (ramp_init_y_up) {
			ramps[4] = ramp1d(y, y_rampup_start, y_rampup_end, false);
		}
		if (ramp_init_y_down) {
			ramps[5] = ramp1d(y, y_rampdown_start, y_rampdown_end, true);
		}
		int min = *std::min_element(std::begin(ramps), std::end(ramps));
		return ramps[min];
	}
	else {
		return 1.0;
	}
}