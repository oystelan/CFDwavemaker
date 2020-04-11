#include "sgrid.h"
#include "omp.h"
#include <algorithm>
#include <iostream>
#include <cmath>

#define PI 3.1415926535897


// A streching function for setting variable layer thickness
double sGrid::slayer(int layerno) {
	//fprintf(stdout,"numlayers: %d",nl);
	double sfac = 3.0; // fixme: this should be made dimensionless and a function of specified wave
	double* dd = new double[nl];
	double ddsum = 0.;
	for (int ii = 0; ii < nl; ii++) {
		//fprintf(stdout,"%d",ii);
		dd[ii] = (1. + (nl - ii) * sfac);
		ddsum += dd[ii];
	}
	double layer_percentage = dd[layerno] / ddsum;
	delete[] dd;
	return layer_percentage;
}

// A function for constant layer thickness (equal spacing as a function of z)
double sGrid::clayer(int layerno) {
	//fprintf(stdout,"numlayers: %d",nl);
	return 1./double(nl);
}

// Transforms normal z axis to streched sigma coordinates 
double sGrid::z2s(double z, double wave_elev, double depth) {
	return (z - wave_elev) / (depth + wave_elev);
}

// Transforms stretched sigma coordinate to normal z
double sGrid::s2z(double s, double wave_elev, double depth) {
	return wave_elev + s * (depth + wave_elev);
}



double sGrid::s2tan(double s) {
	// s defined between 0 and 1, where 0 is seabed, 1 is sea surface
	// returns tangens strethced coordintates tan, which is also defined between 0 and 1
	double a_end = 70;
	double a_potens = 1.5;

	return 1. - (std::pow(std::tan(s * a_end * PI / 180.), a_potens) / std::pow(std::tan(a_end * PI / 180.), a_potens));

}

double sGrid::tan2s(double t) {
	// The inverse of the above function s2tan. from tan stretched to normal constant spacing
	double a_end = 70;
	double a_potens = 1.5;
	return std::atan(std::pow((1. - t) * std::pow(std::tan(a_end * PI / 180.), a_potens), (1 / a_potens))) / (a_end *PI / 180.);

}

// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 3
void sGrid::initialize_kinematics(Irregular& irregular, double tpt) {
	// Allocating memory for storage of surface elevation and velocities
	UX0 = new double[NX * NY * NZ];
	UY0 = new double[NX * NY * NZ];
	UZ0 = new double[NX * NY * NZ];

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
				eta_temp = irregular.eta(tpt, xpt, ypt);

				double Ux0 = irregular.u1(tpt, xpt, ypt, 0.0) + irregular.u2(tpt, xpt, ypt, 0.0);
				double Uy0 = irregular.v1(tpt, xpt, ypt, 0.0) + irregular.v2(tpt, xpt, ypt, 0.0);
				double Uz0 = irregular.w1(tpt, xpt, ypt, 0.0) + irregular.w2(tpt, xpt, ypt, 0.0);

				double PHI_dxdz = irregular.phi1_dxdz(tpt, xpt, ypt);
				double PHI_dydz = irregular.phi1_dydz(tpt, xpt, ypt);
				double PHI_dzdz = irregular.phi1_dzdz(tpt, xpt, ypt);

				for (int m = 0; m < NZ; m++) {
					zpt = domain_start[2] + dz * m;
					if (zpt > (eta_temp + dz)) {
						UX0[i * NY * NZ + j * NZ + m] = 0.0;
						UY0[i * NY * NZ + j * NZ + m] = 0.0;
						UZ0[i * NY * NZ + j * NZ + m] = 0.0;
					}
					else if (zpt > 0.) {
						UX0[i * NY * NZ + j * NZ + m] = Ux0 + PHI_dxdz * zpt;
						UY0[i * NY * NZ + j * NZ + m] = Uy0 + PHI_dydz * zpt;
						UZ0[i * NY * NZ + j * NZ + m] = Uz0 + PHI_dzdz * zpt;
					}
					else {
						UX0[i * NY * NZ + j * NZ + m] = irregular.u1(tpt, xpt, ypt, zpt) + irregular.u2(tpt, xpt, ypt, zpt);
						UY0[i * NY * NZ + j * NZ + m] = irregular.v1(tpt, xpt, ypt, zpt) + irregular.v2(tpt, xpt, ypt, zpt);
						UZ0[i * NY * NZ + j * NZ + m] = irregular.w1(tpt, xpt, ypt, zpt) + irregular.w2(tpt, xpt, ypt, zpt);
					}
					/*UX[i*NY*NZ + j*NZ + m] = uu(tpt, xpt, ypt, zpt);
					UY[i*NY*NZ + j*NZ + m] = vv(tpt, xpt, ypt, zpt);
					UZ[i*NY*NZ + j*NZ + m] = ww(tpt, xpt, ypt, zpt);*/
				}
			}
		}
	} // End parallel initialization

	std::cout << "Generation of domain kinematics data completed. ";
	dd = omp_get_wtime() - dd;
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
	std::cout << "Interpolation can commence..." << std::endl;
	initkin = 1;
}


// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 3
void sGrid::initialize_surface_elevation(Irregular& irregular, double tpt) {

	// Allocating memory for storage of surface elevation and velocities
	ETA0 = new double[NX * NY];

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
				ETA0[i * NY + j] = irregular.eta1(tpt, xpt, ypt) + irregular.eta2(tpt, xpt, ypt);
			}
		}
	}
	dd = omp_get_wtime() - dd;

	std::cout << "Surface Elevation generated successfully. ";
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
	initsurf = 1;
}

/* Function for trilinear interpolation on a cartesian evenly spaced mesh*/
double sGrid::trilinear_interpolation(double* VAR, double xpt, double ypt, double zpt, int _nx, int _ny, int _nz, double _dx, double _dy, double _dz, double* domain) {
	double nxp = std::min(double(_nx), std::max(0., (xpt - domain[0]) / _dx));
	double nyp = std::min(double(_ny), std::max(0., (ypt - domain[1]) / _dy));
	double nzp = std::min(double(_nz), std::max(0., (zpt - domain[2]) / _dz));

	double C000 = VAR[int(floor(nxp) * _ny * _nz + floor(nyp) * _nz + floor(nzp))];
	double C001 = VAR[int(floor(nxp) * _ny * _nz + floor(nyp) * _nz + ceil(nzp))];
	double C010 = VAR[int(floor(nxp) * _ny * _nz + ceil(nyp) * _nz + floor(nzp))];
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

/* bilinear interpolation function used to interpolate surface values on a regular evenly spaced grid*/
double sGrid::bilinear_interpolation(double* VAR, double xpt, double ypt, int _nx, int _ny, double _dx, double _dy, double* domain) {

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


double sGrid::u(double xpt, double ypt, double zpt) {
	return trilinear_interpolation(UX0, xpt, ypt, z2s(zpt, eta(xpt, ypt), water_depth), NX, NY, NZ, dx, dy, dz, domain_start);
}

double sGrid::v(double xpt, double ypt, double zpt) {
	return trilinear_interpolation(UY0, xpt, ypt, z2s(zpt, eta(xpt, ypt), water_depth), NX, NY, NZ, dx, dy, dz, domain_start);
}

double sGrid::w(double xpt, double ypt, double zpt) {
	return trilinear_interpolation(UZ0, xpt, ypt, z2s(zpt, eta(xpt, ypt), water_depth), NX, NY, NZ, dx, dy, dz, domain_start);
}

double sGrid::eta(double xpt, double ypt) {
	return bilinear_interpolation(ETA0, xpt, ypt, NX, NY, dx, dy, domain_start);
}


bool sGrid::CheckTime(double tpt) {
	/* Checks to see if the time tpt is within the interval t0 to t1. If so, returns true*/
	if (tpt > t1) {
		return false;
		std::cout << "T0: " << t0 << ", t1: " << t1 << ", tpt: " << tpt << std::endl;
	}
	else {
		return true;
	}
}

// function to find if given point 
// lies inside a given rectangle or not. 
bool sGrid::CheckBounds(double* bounds, double x, double y, double z)
{
	if (x >= bounds[0] && x <= bounds[1] && y >= bounds[2] && y <= bounds[3] && z >= bounds[4] && z <= bounds[5])
		return true;
	else {
		std::cout << "position x,y,z=" << x << "," << y << "," << z << " out of bounds. updating wallx borders" << std::endl;
		return false;
	}
}