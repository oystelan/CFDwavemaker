#include "sgrid.h"
#include "omp.h"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <string>
#if defined(_WIN32)
#include < direct.h >
#endif


//----------------------------------------------------------------------------------------------------------------------------------------


// LSgrid functions


//----------------------------------------------------------------------------------------------------------------------------------------



// Allocation of memory to storage matrices
void sGrid::allocate() {
	ETA0 = new double[nx * ny];
	ETA1 = new double[nx * ny];

	IGNORE = new int[nx * ny];

	UX0 = new double[nx * ny * nl];
	UY0 = new double[nx * ny * nl];
	UZ0 = new double[nx * ny * nl];

	UX1 = new double[nx * ny * nl];
	UY1 = new double[nx * ny * nl];
	UZ1 = new double[nx * ny * nl];
}

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
// s defined between -1 (seabed) and 0 (free surface)
double sGrid::z2s(double z, double wave_elev, double depth) {
	return (z - wave_elev) / (depth + wave_elev);
}

// Transforms stretched sigma coordinate to normal z
double sGrid::s2z(double s, double wave_elev, double depth) {
	return wave_elev + s * (depth + wave_elev);
}

double sGrid::s2tan(double s) {
	// s defined between -1 and 0, where -1 is seabed, 0 is sea surface
	// returns tangens strethced coordintates tan, which is also defined between -1 and 0	
	return -std::pow(std::tan(-s * tan_a) , tan_b) / std::pow(std::tan(tan_a), tan_b);
}

double sGrid::tan2s(double t) {
	// The inverse of the above function s2tan. from tan stretched to normal constant spacing	
	return -std::atan(std::pow(-t * std::pow(std::tan(tan_a) , tan_b) , 1. / tan_b)) / tan_a;
}

/* Function for trilinear interpolation on a cartesian evenly spaced mesh*/
double sGrid::trilinear_interpolation(double* VAR0, double* VAR1, double tpt, double xpt, double ypt, double zpt) {

	double nxp = std::min(double(nx), std::max(0., (xpt - domain[0]) / dx));
	double nyp = std::min(double(ny), std::max(0., (ypt - domain[2]) / dy));

	// Find surface elevation (required for stretching)
	double C00 = ETA0[int(floor(nxp) * ny + floor(nyp))];
	double C01 = ETA0[int(floor(nxp) * ny + ceil(nyp))];
	double C10 = ETA0[int(ceil(nxp) * ny + floor(nyp))];
	double C11 = ETA0[int(ceil(nxp) * ny + ceil(nyp))];
	double D00 = ETA1[int(floor(nxp) * ny + floor(nyp))];
	double D01 = ETA1[int(floor(nxp) * ny + ceil(nyp))];
	double D10 = ETA1[int(ceil(nxp) * ny + floor(nyp))];
	double D11 = ETA1[int(ceil(nxp) * ny + ceil(nyp))];

	double xd = nxp - floor(nxp);
	double yd = nyp - floor(nyp);

	double td = std::min(1., std::max(0., (tpt - t0) / dt));

	double C0 = C00 * (1. - xd) + C10 * xd;
	double C1 = C01 * (1. - xd) + C11 * xd;
	double D0 = D00 * (1. - xd) + D10 * xd;
	double D1 = D01 * (1. - xd) + D11 * xd;

	double wave_elev0 = C0 * (1. - yd) + C1 * yd;
	double wave_elev1 = D0 * (1. - yd) + D1 * yd;

	/* Todo: fix this check at some point
	if (zpt > wave_elev) {
		return 0.;
	}*/

	//std::cout << "dy: " << dy << ", nyp: " << nyp << ", ypt: " << ypt << std::endl;

	double spt0 = tan2s(std::max(z2s(std::min(zpt - swl, wave_elev0 - swl), wave_elev0 - swl, water_depth), -1.));
	double nsp0 = (spt0 + 1.) / ds;
	double spt1 = tan2s(std::max(z2s(std::min(zpt - swl, wave_elev1 - swl), wave_elev1 - swl, water_depth), -1.));
	double nsp1 = (spt1 + 1.) / ds;

	//std::cout << "orig: nsp0:" << nsp0 << ", nsp1: " << nsp1 << "spt0: " << spt0 << ", wavelev0: " << wave_elev0 << std::endl;
	//exit(0);
	// Trilinear interpolation.
	double C000 = VAR0[int(floor(nxp) * ny * nl + floor(nyp) * nl + floor(nsp0))];
	double C001 = VAR0[int(floor(nxp) * ny * nl + floor(nyp) * nl + ceil(nsp0))];
	double C010 = VAR0[int(floor(nxp) * ny * nl + ceil(nyp) * nl + floor(nsp0))];
	double C011 = VAR0[int(floor(nxp) * ny * nl + ceil(nyp) * nl + ceil(nsp0))];
	double C100 = VAR0[int(ceil(nxp) * ny * nl + floor(nyp) * nl + floor(nsp0))];
	double C101 = VAR0[int(ceil(nxp) * ny * nl + floor(nyp) * nl + ceil(nsp0))];
	double C110 = VAR0[int(ceil(nxp) * ny * nl + ceil(nyp) * nl + floor(nsp0))];
	double C111 = VAR0[int(ceil(nxp) * ny * nl + ceil(nyp) * nl + ceil(nsp0))];
	double D000 = VAR1[int(floor(nxp) * ny * nl + floor(nyp) * nl + floor(nsp1))];
	double D001 = VAR1[int(floor(nxp) * ny * nl + floor(nyp) * nl + ceil(nsp1))];
	double D010 = VAR1[int(floor(nxp) * ny * nl + ceil(nyp) * nl + floor(nsp1))];
	double D011 = VAR1[int(floor(nxp) * ny * nl + ceil(nyp) * nl + ceil(nsp1))];
	double D100 = VAR1[int(ceil(nxp) * ny * nl + floor(nyp) * nl + floor(nsp1))];
	double D101 = VAR1[int(ceil(nxp) * ny * nl + floor(nyp) * nl + ceil(nsp1))];
	double D110 = VAR1[int(ceil(nxp) * ny * nl + ceil(nyp) * nl + floor(nsp1))];
	double D111 = VAR1[int(ceil(nxp) * ny * nl + ceil(nyp) * nl + ceil(nsp1))];

	double sd0 = nsp0 - floor(nsp0);
	double sd1 = nsp1 - floor(nsp1);

	C00 = C000 * (1. - xd) + C100 * xd;
	C01 = C001 * (1. - xd) + C101 * xd;
	C10 = C010 * (1. - xd) + C110 * xd;
	C11 = C011 * (1. - xd) + C111 * xd;
	//std::cout << int(ceil(nxp) * ny * nl + ceil(nyp) * nl + ceil(nsp0)) << ", " << ceil(nxp) << ", " << ceil(nyp) << ", " << ceil(nsp0) << std::endl;

	//std::cout << C010 << ", " << C110 << ", " << C011 << ", " << C111 << std::endl;

	D00 = D000 * (1. - xd) + D100 * xd;
	D01 = D001 * (1. - xd) + D101 * xd;
	D10 = D010 * (1. - xd) + D110 * xd;
	D11 = D011 * (1. - xd) + D111 * xd;

	C0 = C00 * (1. - yd) + C10 * yd;
	C1 = C01 * (1. - yd) + C11 * yd;
	D0 = D00 * (1. - yd) + D10 * yd;
	D1 = D01 * (1. - yd) + D11 * yd;

	return (C0 * (1. - sd0) + C1 * sd0) * (1. - td) + (D0 * (1. - sd1) + D1 * sd1) * td;
}
template <typename T>
T clip(const T& n, const T& lower, const T& upper) {
	return std::max(lower, std::min(n, upper));
}

/* Function for trilinear interpolation on a cartesian evenly spaced mesh*/
double sGrid::trilinear_interpolation2(double* VAR0, double* VAR1, double tpt, double xpt, double ypt, double zpt) {

	float nxp_temp, nyp_temp;
	double xd = std::modf(clip((xpt - domain[0]) / dx, 0., nx - 1.), &nxp_temp);
	double yd = std::modf(clip((ypt - domain[2]) / dy, 0., ny - 1.), &nyp_temp);

	int nxp = int(nxp_temp);
	int nyp = int(nyp_temp);

	double C00 = ETA0[nxp * ny + nyp];
	double C01 = ETA0[nxp * ny + clip(nyp + 1, 0, ny - 1)];
	double C10 = ETA0[clip(nxp + 1, 0, nx - 1) * ny + nyp];
	double C11 = ETA0[clip(nxp + 1, 0, nx - 1) * ny + clip(nyp + 1, 0, ny - 1)];
	double D00 = ETA1[nxp * ny + nyp];
	double D01 = ETA1[nxp * ny + clip(nyp + 1, 0, ny - 1)];
	double D10 = ETA1[clip(nxp + 1, 0, nx - 1) * ny + nyp];
	double D11 = ETA1[clip(nxp + 1, 0, nx - 1) * ny + clip(nyp + 1, 0, ny - 1)];

	double td = std::min(1., std::max(0., (tpt - t0) / dt));

	double C0 = C00 * (1. - xd) + C10 * xd;
	double C1 = C01 * (1. - xd) + C11 * xd;
	double D0 = D00 * (1. - xd) + D10 * xd;
	double D1 = D01 * (1. - xd) + D11 * xd;

	double wave_elev0 = C0 * (1. - yd) + C1 * yd;
	double wave_elev1 = D0 * (1. - yd) + D1 * yd;

	/* Todo: fix this check at some point
	if (zpt > wave_elev) {
		return 0.;
	}*/

	//std::cout << "dy: " << dy << ", nyp: " << nyp << ", ypt: " << ypt << std::endl;

	double spt0 = tan2s(std::max(z2s(std::min(zpt - swl, wave_elev0 - swl), wave_elev0 - swl, water_depth), -1.));
	double nsp0a = (spt0 + 1.) / ds;
	double spt1 = tan2s(std::max(z2s(std::min(zpt - swl, wave_elev1 - swl), wave_elev1 - swl, water_depth), -1.));
	double nsp1a = (spt1 + 1.) / ds;

	float nsp0_temp, nsp1_temp;
	double sd0 = std::modf((spt0 + 1.) / ds, &nsp0_temp);
	double sd1 = std::modf((spt1 + 1.) / ds, &nsp1_temp);
	int nsp0 = int(nsp0_temp);
	int nsp1 = int(nsp1_temp);

	//std::cout << floor(nsp0a) << " " << nsp0 << "spt0: " << spt0 << ", wavelev0: " << wave_elev0 << std::endl;

	//std::cout << "nsp0:" << nsp0 << ", nsp1: " << nsp1 << std::endl;

	//exit(0);
	// Trilinear interpolation.
	double C000 = VAR0[nxp * ny * nl + nyp * nl + nsp0];
	double C001 = VAR0[nxp * ny * nl + nyp * nl + clip(nsp0 + 1, 0, nl - 1)];
	double C010 = VAR0[nxp * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + nsp0];
	double C011 = VAR0[nxp * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + clip(nsp0 + 1, 0, nl - 1)];
	double C100 = VAR0[clip(nxp + 1, 0, nx - 1) * ny * nl + nyp * nl + nsp0];
	double C101 = VAR0[clip(nxp + 1, 0, nx - 1) * ny * nl + nyp * nl + clip(nsp0 + 1, 0, nl - 1)];
	double C110 = VAR0[clip(nxp + 1, 0, nx - 1) * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + nsp0];//
	double C111 = VAR0[clip(nxp + 1, 0, nx - 1) * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + clip(nsp0 + 1, 0, nl - 1)];//
	double D000 = VAR1[nxp * ny * nl + nyp * nl + nsp1];
	double D001 = VAR1[nxp * ny * nl + nyp * nl + clip(nsp1 + 1, 0, nl - 1)];
	double D010 = VAR1[nxp * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + nsp1];
	double D011 = VAR1[nxp * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + clip(nsp1 + 1, 0, nl - 1)];
	double D100 = VAR1[clip(nxp + 1, 0, nx - 1) * ny * nl + nyp * nl + nsp1];
	double D101 = VAR1[clip(nxp + 1, 0, nx - 1) * ny * nl + nyp * nl + clip(nsp1 + 1, 0, nl - 1)];
	double D110 = VAR1[clip(nxp + 1, 0, nx - 1) * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + nsp1];
	double D111 = VAR1[clip(nxp + 1, 0, nx - 1) * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + clip(nsp1 + 1, 0, nl - 1)];

	//double sd0 = nsp0 - floor(nsp0);
	//double sd1 = nsp1 - floor(nsp1);

	C00 = C000 * (1. - xd) + C100 * xd;
	C01 = C001 * (1. - xd) + C101 * xd;
	C10 = C010 * (1. - xd) + C110 * xd;
	C11 = C011 * (1. - xd) + C111 * xd;
	//std::cout << int(ceil(nxp) * ny * nl + ceil(nyp) * nl + ceil(nsp0)) << ", " << ceil(nxp) << ", " << ceil(nyp) << ", " << ceil(nsp0) << std::endl;

	//std::cout << C010 << ", " << C110 << ", " << C011 << ", " << C111 << std::endl;

	D00 = D000 * (1. - xd) + D100 * xd;
	D01 = D001 * (1. - xd) + D101 * xd;
	D10 = D010 * (1. - xd) + D110 * xd;
	D11 = D011 * (1. - xd) + D111 * xd;

	C0 = C00 * (1. - yd) + C10 * yd;
	C1 = C01 * (1. - yd) + C11 * yd;
	D0 = D00 * (1. - yd) + D10 * yd;
	D1 = D01 * (1. - yd) + D11 * yd;

	return (C0 * (1. - sd0) + C1 * sd0) * (1. - td) + (D0 * (1. - sd1) + D1 * sd1) * td;
}

/* bilinear interpolation function used to interpolate surface values on a regular evenly spaced grid*/
double sGrid::bilinear_interpolation(double* VAR0, double* VAR1, double tpt, double xpt, double ypt) {

	double nxp = std::min(double(nx), std::max(0., (xpt - domain[0]) / dx));
	double nyp = std::min(double(ny), std::max(0., (ypt - domain[2]) / dy));


	double C00 = VAR0[int(floor(nxp) * ny + floor(nyp))];
	double C01 = VAR0[int(floor(nxp) * ny + ceil(nyp))];
	double C10 = VAR0[int(ceil(nxp) * ny + floor(nyp))];
	double C11 = VAR0[int(ceil(nxp) * ny + ceil(nyp))];
	double D00 = VAR1[int(floor(nxp) * ny + floor(nyp))];
	double D01 = VAR1[int(floor(nxp) * ny + ceil(nyp))];
	double D10 = VAR1[int(ceil(nxp) * ny + floor(nyp))];
	double D11 = VAR1[int(ceil(nxp) * ny + ceil(nyp))];

	double xd = nxp - floor(nxp);
	double yd = nyp - floor(nyp);
	double td = std::min(1., std::max(0., (tpt - t0) / dt));

	double C0 = C00 * (1. - xd) + C10 * xd;
	double C1 = C01 * (1. - xd) + C11 * xd;
	double D0 = D00 * (1. - xd) + D10 * xd;
	double D1 = D01 * (1. - xd) + D11 * xd;

	return (C0 * (1. - yd) + C1 * yd) * (1. - td) + (D0 * (1. - yd) + D1 * yd) * td;
}

//#define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))



/* bilinear interpolation function used to interpolate surface values on a regular evenly spaced grid*/
double sGrid::bilinear_interpolation2(double* VAR0, double* VAR1, double tpt, double xpt, double ypt) {



	float nxp_temp, nyp_temp;
	double xd = std::modf(clip((xpt - domain[0]) / dx, 0., nx - 1.), &nxp_temp);
	double yd = std::modf(clip((ypt - domain[2]) / dy, 0., ny - 1.), &nyp_temp);

	int nxp = int(nxp_temp);
	int nyp = int(nyp_temp);

	//std::cout << nxp << " " << nyp << "VAR0: " << VAR0[0] << std::endl;

	double C00 = VAR0[nxp * ny + nyp];
	double C01 = VAR0[nxp * ny + clip(nyp + 1, 0, ny - 1)];
	double C10 = VAR0[clip(nxp + 1, 0, nx - 1) * ny + nyp];
	double C11 = VAR0[clip(nxp + 1, 0, nx - 1) * ny + clip(nyp + 1, 0, ny - 1)];
	double D00 = VAR1[nxp * ny + nyp];
	double D01 = VAR1[nxp * ny + clip(nyp + 1, 0, ny - 1)];
	double D10 = VAR1[clip(nxp + 1, 0, nx - 1) * ny + nyp];
	double D11 = VAR1[clip(nxp + 1, 0, nx - 1) * ny + clip(nyp + 1, 0, ny - 1)];

	//std::cout << C00 << " " << C01 << " " << C10 << " " << C11 << " " << D00 << " " << D01 << " " << D10 << " " << D11 << std::endl;
	//double xd = nxp - floor(nxp);
	//double yd = nyp - floor(nyp);
	double td = std::min(1., std::max(0., (tpt - t0) / dt));

	double C0 = C00 * (1. - xd) + C10 * xd;
	double C1 = C01 * (1. - xd) + C11 * xd;
	double D0 = D00 * (1. - xd) + D10 * xd;
	double D1 = D01 * (1. - xd) + D11 * xd;

	return (C0 * (1. - yd) + C1 * yd) * (1. - td) + (D0 * (1. - yd) + D1 * yd) * td;
}


double sGrid::u(double tpt, double xpt, double ypt, double zpt) {
	return trilinear_interpolation2(UX0, UX1, tpt, xpt, ypt, zpt);
}

double sGrid::v(double tpt, double xpt, double ypt, double zpt) {
	return trilinear_interpolation2(UY0, UY1, tpt, xpt, ypt, zpt);
}

double sGrid::w(double tpt, double xpt, double ypt, double zpt) {
	return trilinear_interpolation2(UZ0, UZ1, tpt, xpt, ypt, zpt);
}

double sGrid::eta(double tpt, double xpt, double ypt) {
	//update_bounds(xpt, ypt);
	return bilinear_interpolation2(ETA0, ETA1, tpt, xpt, ypt);
}


bool sGrid::CheckTime(double tpt) {
	/* Checks to see if the time tpt is within the interval t0 to t1. If so, returns true*/
	if (tpt > t0 + dt) {
		std::cout << "t0: " << t0 << ", t1: " << (t0 + dt) << ", tpt: " << tpt << std::endl;
		return false;

	}
	return true;
}

// function to find if given point 
// lies inside a given rectangle or not. 
bool sGrid::CheckBounds()
{
	if (bxmin >= domain[0] && bxmax <= domain[1] && bymin >= domain[2] && bymax <= domain[3])
		return true;
	else {
		std::cout << "Requested point outside specified grid domain. adjust the bounds of the grid and try again." << std::endl;
		exit(-1);
		return false;
	}
}

void sGrid::update_bounds(double xpt, double ypt) {
	bxmin = std::min(xpt, bxmin);
	bxmax = std::min(xpt, bxmax);
	bymin = std::min(ypt, bymin);
	bymax = std::min(ypt, bymax);

}


//----------------------------------------------------------------------------------------------------------------------------------------


// Second order theory functions


//----------------------------------------------------------------------------------------------------------------------------------------

// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 3
void sGrid::initialize_kinematics(Irregular& irregular) {

	dx = (domain[1] - domain[0]) / std::max(1., double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));
	//dz = (domain_end[2] - domain_start[2]) / double(NZ - 1);
	ds = 1. / std::max(1., double(nl - 1));

	double dd = omp_get_wtime();

	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());

#pragma omp parallel // start parallel initialization
	{
#pragma omp master
		std::cout << "Number of available threads: " << omp_get_num_threads() << std::endl;

		// Main grid
#pragma omp for collapse(2)
		for (int i = 0; i < nx; i++) {

			for (int j = 0; j < ny; j++) {
				double xpt = domain[0] + dx * i;
				double ypt = domain[2] + dy * j;
				double eta0_temp = ETA0[i * ny + j] - swl;

				double PHI0_dxdz = irregular.phi1_dxdz(t0, xpt, ypt);
				double PHI0_dydz = irregular.phi1_dydz(t0, xpt, ypt);
				double PHI0_dzdz = irregular.phi1_dzdz(t0, xpt, ypt);

				for (int m = 0; m < nl; m++) {
					double spt = s2tan(-1. + ds * m);
					double zpt0 = s2z(spt, eta0_temp, water_depth);

					UX0[i * ny * nl + j * nl + m] = irregular.u1(t0, xpt, ypt, zpt0) + irregular.u2(t0, xpt, ypt, zpt0) + PHI0_dxdz * std::max(0., (zpt0 - swl));
					UY0[i * ny * nl + j * nl + m] = irregular.v1(t0, xpt, ypt, zpt0) + irregular.v2(t0, xpt, ypt, zpt0) + PHI0_dydz * std::max(0., (zpt0 - swl));
					UZ0[i * ny * nl + j * nl + m] = irregular.w1(t0, xpt, ypt, zpt0) + irregular.w2(t0, xpt, ypt, zpt0) + PHI0_dzdz * std::max(0., (zpt0 - swl));
				}
			}
		}
	}
		
	if (init_only) {
		UX1 = UX0;
		UY1 = UY0;
		UZ1 = UZ0;
	}
	else{
	#pragma omp parallel // start parallel initialization
		{
#pragma omp for collapse(2)
			for (int i = 0; i < nx; i++) {
				for (int j = 0; j < ny; j++) {
					double xpt = domain[0] + dx * i;
					double ypt = domain[2] + dy * j;
					double eta1_temp = ETA1[i * ny + j] - swl;

					double PHI1_dxdz = irregular.phi1_dxdz(t0 + dt, xpt, ypt);
					double PHI1_dydz = irregular.phi1_dydz(t0 + dt, xpt, ypt);
					double PHI1_dzdz = irregular.phi1_dzdz(t0 + dt, xpt, ypt);

					for (int m = 0; m < nl; m++) {
						double spt = s2tan(-1. + ds * m);
						double zpt1 = s2z(spt, eta1_temp, water_depth);

						UX1[i * ny * nl + j * nl + m] = irregular.u1(t0 + dt, xpt, ypt, zpt1) + irregular.u2(t0 + dt, xpt, ypt, zpt1) + PHI1_dxdz * std::max(0., (zpt1 - swl));
						UY1[i * ny * nl + j * nl + m] = irregular.v1(t0 + dt, xpt, ypt, zpt1) + irregular.v2(t0 + dt, xpt, ypt, zpt1) + PHI1_dydz * std::max(0., (zpt1 - swl));
						UZ1[i * ny * nl + j * nl + m] = irregular.w1(t0 + dt, xpt, ypt, zpt1) + irregular.w2(t0 + dt, xpt, ypt, zpt1) + PHI1_dzdz * std::max(0., (zpt1 - swl));
					}
				}
			}
		}
	} 
	// End parallel initialization
	if (dump_vtk) {
		write_vtk(false);
		update_count++;
	}
	std::cout << "Generation of domain kinematics data using irregular second order wave theory completed. ";
	dd = omp_get_wtime() - dd;
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
	std::cout << "Interpolation can commence..." << std::endl;
}

// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 3
void sGrid::initialize_kinematics_with_ignore(Irregular& irregular) {

	dx = (domain[1] - domain[0]) / std::max(1., double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));
	//dz = (domain_end[2] - domain_start[2]) / double(NZ - 1);
	ds = 1. / std::max(1., double(nl - 1));

	double dd = omp_get_wtime();

	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());

#pragma omp parallel // start parallel initialization
	{
#pragma omp master
		std::cout << "Number of available threads: " << omp_get_num_threads() << std::endl;

		// Main grid
#pragma omp for collapse(2)
		for (int i = 0; i < nx; i++) {

			for (int j = 0; j < ny; j++) {
				if (!IGNORE[i * ny + j]) {
					double xpt = domain[0] + dx * i;
					double ypt = domain[2] + dy * j;
					double eta0_temp = ETA0[i * ny + j] - swl;

					double PHI0_dxdz = irregular.phi1_dxdz(t0, xpt, ypt);
					double PHI0_dydz = irregular.phi1_dydz(t0, xpt, ypt);
					double PHI0_dzdz = irregular.phi1_dzdz(t0, xpt, ypt);

					for (int m = 0; m < nl; m++) {
						double spt = s2tan(-1. + ds * m);
						double zpt0 = s2z(spt, eta0_temp, water_depth);

						UX0[i * ny * nl + j * nl + m] = irregular.u1(t0, xpt, ypt, zpt0) + irregular.u2(t0, xpt, ypt, zpt0) + PHI0_dxdz * std::max(0., (zpt0 - swl));
						UY0[i * ny * nl + j * nl + m] = irregular.v1(t0, xpt, ypt, zpt0) + irregular.v2(t0, xpt, ypt, zpt0) + PHI0_dydz * std::max(0., (zpt0 - swl));
						UZ0[i * ny * nl + j * nl + m] = irregular.w1(t0, xpt, ypt, zpt0) + irregular.w2(t0, xpt, ypt, zpt0) + PHI0_dzdz * std::max(0., (zpt0 - swl));
					}
				}
			}
		}
	}

	if (init_only) {
		UX1 = UX0;
		UY1 = UY0;
		UZ1 = UZ0;
	}
	else {
#pragma omp parallel // start parallel initialization
		{
#pragma omp for collapse(2)
			for (int i = 0; i < nx; i++) {
				for (int j = 0; j < ny; j++) {
					if (!IGNORE[i * ny + j]) {
						double xpt = domain[0] + dx * i;
						double ypt = domain[2] + dy * j;
						double eta1_temp = ETA1[i * ny + j] - swl;

						double PHI1_dxdz = irregular.phi1_dxdz(t0 + dt, xpt, ypt);
						double PHI1_dydz = irregular.phi1_dydz(t0 + dt, xpt, ypt);
						double PHI1_dzdz = irregular.phi1_dzdz(t0 + dt, xpt, ypt);

						for (int m = 0; m < nl; m++) {
							double spt = s2tan(-1. + ds * m);
							double zpt1 = s2z(spt, eta1_temp, water_depth);

							UX1[i * ny * nl + j * nl + m] = irregular.u1(t0 + dt, xpt, ypt, zpt1) + irregular.u2(t0 + dt, xpt, ypt, zpt1) + PHI1_dxdz * std::max(0., (zpt1 - swl));
							UY1[i * ny * nl + j * nl + m] = irregular.v1(t0 + dt, xpt, ypt, zpt1) + irregular.v2(t0 + dt, xpt, ypt, zpt1) + PHI1_dydz * std::max(0., (zpt1 - swl));
							UZ1[i * ny * nl + j * nl + m] = irregular.w1(t0 + dt, xpt, ypt, zpt1) + irregular.w2(t0 + dt, xpt, ypt, zpt1) + PHI1_dzdz * std::max(0., (zpt1 - swl));
						}
					}
				}
			}
		}
	}
	// End parallel initialization
	if (dump_vtk) {
		write_vtk(false);
		update_count++;
	}
	std::cout << "Generation of domain kinematics data using irregular second order wave theory completed. ";
	dd = omp_get_wtime() - dd;
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
	std::cout << "Interpolation can commence..." << std::endl;
}



// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 4
void sGrid::initialize_surface_elevation(Irregular& irregular, double t_target) {
	std::cout << "time: " << t_target << std::endl;
	t0 = t_target;

	// Allocating memory for storage of surface elevation and velocities

	dx = (domain[1] - domain[0]) / std::max(1., double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));

	bxmin = domain[0];
	bxmax = domain[1];
	bymin = domain[2];
	bymax = domain[3];

	double dd = omp_get_wtime();
	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());

#pragma omp parallel
	{
		// Main grid
#pragma omp for collapse(2)
		for (int i = 0; i < nx; i++) {

			for (int j = 0; j < ny; j++) {
				double xpt = domain[0] + dx * i;
				double ypt = domain[2] + dy * j;

				ETA0[i * ny + j] = irregular.eta1(t0, xpt, ypt) + irregular.eta2(t0, xpt, ypt) + swl;
				ETA1[i * ny + j] = irregular.eta1(t0 + dt, xpt, ypt) + irregular.eta2(t0 + dt, xpt, ypt) + swl;

			}
		}
	}
	dd = omp_get_wtime() - dd;

	std::cout << "Surface Elevation generated successfully. ";
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
}
// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 3
void sGrid::initialize_surface_elevation_with_ignore(Irregular& irregular, double t_target) {

	std::cout << "time: " << t_target << std::endl;
	t0 = t_target;

	// Allocating memory for storage of surface elevation and velocities

	dx = (domain[1] - domain[0]) / std::max(1., double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));

	bxmin = domain[0];
	bxmax = domain[1];
	bymin = domain[2];
	bymax = domain[3];

	double dd = omp_get_wtime();
	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());

#pragma omp parallel
	{
		// Main grid
#pragma omp for collapse(2)
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				double xpt = domain[0] + dx * i;
				double ypt = domain[2] + dy * j;
				if (!IGNORE[i * ny + j]) {
					ETA0[i * ny + j] = irregular.eta1(t0, xpt, ypt) + irregular.eta2(t0, xpt, ypt) + swl;
					ETA1[i * ny + j] = irregular.eta1(t0 + dt, xpt, ypt) + irregular.eta2(t0 + dt, xpt, ypt) + swl;
				}
				else {
					ETA0[i * ny + j] = swl;
					ETA1[i * ny + j] = swl;
				}
			}
		}
	}
	dd = omp_get_wtime() - dd;

	std::cout << "Surface Elevation generated successfully. ";
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
}


// When called, updates the arrays storing surface elevation and kinematics data for timestep t0 = t1, t1 = t1+dt
void sGrid::update(Irregular& irregular, double t_target)
{
	// Start by checking bounds
	/*
	if (!disable_checkbounds){
		CheckBounds();
	}*/

	// new time step
	if ((t_target / dt - (t0 + 2 * dt) / dt) > 0.) {
		double new_time = dt * std::floor(t_target / dt);
		initialize_surface_elevation(irregular, new_time);
		initialize_kinematics(irregular);
	}
	else {
		t0 += dt;
		// Updating surface elevations
		double dd = omp_get_wtime();
		//omp_set_num_threads(1);
		omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel
		{
			// Main grid
#pragma omp for collapse(2)
			for (int i = 0; i < nx; i++) {
				for (int j = 0; j < ny; j++) {
					double xpt = domain[0] + dx * i;
					//std::cout << "processornum: " << omp_get_thread_num() << std::endl;
					double ypt = domain[2] + dy * j;
					if (!IGNORE[i * ny + j]) {
						ETA0[i * ny + j] = ETA1[i * ny + j];
						ETA1[i * ny + j] = irregular.eta1(t0 + dt, xpt, ypt) + irregular.eta2(t0 + dt, xpt, ypt) + swl;

						double Ux1 = irregular.u1(t0 + dt, xpt, ypt, 0.0) + irregular.u2(t0 + dt, xpt, ypt, 0.0);
						double Uy1 = irregular.v1(t0 + dt, xpt, ypt, 0.0) + irregular.v2(t0 + dt, xpt, ypt, 0.0);
						double Uz1 = irregular.w1(t0 + dt, xpt, ypt, 0.0) + irregular.w2(t0 + dt, xpt, ypt, 0.0);

						double PHI1_dxdz = irregular.phi1_dxdz(t0 + dt, xpt, ypt);
						double PHI1_dydz = irregular.phi1_dydz(t0 + dt, xpt, ypt);
						double PHI1_dzdz = irregular.phi1_dzdz(t0 + dt, xpt, ypt);

						for (int m = 0; m < nl; m++) {
							double spt = s2tan(-1. + ds * m);
							double zpt1 = s2z(spt, ETA1[i * ny + j] - swl, water_depth);

							UX0[i * ny * nl + j * nl + m] = UX1[i * ny * nl + j * nl + m];
							UY0[i * ny * nl + j * nl + m] = UY1[i * ny * nl + j * nl + m];
							UZ0[i * ny * nl + j * nl + m] = UZ1[i * ny * nl + j * nl + m];

							UX1[i * ny * nl + j * nl + m] = irregular.u1(t0 + dt, xpt, ypt, zpt1) + irregular.u2(t0 + dt, xpt, ypt, zpt1) + PHI1_dxdz * std::max(0., (zpt1 - swl));
							UY1[i * ny * nl + j * nl + m] = irregular.v1(t0 + dt, xpt, ypt, zpt1) + irregular.v2(t0 + dt, xpt, ypt, zpt1) + PHI1_dydz * std::max(0., (zpt1 - swl));
							UZ1[i * ny * nl + j * nl + m] = irregular.w1(t0 + dt, xpt, ypt, zpt1) + irregular.w2(t0 + dt, xpt, ypt, zpt1) + PHI1_dzdz * std::max(0., (zpt1 - swl));

						}
					}
				}
			}
		}

		if (dump_vtk) {
			write_vtk(false);
			update_count++;
		}

		dd = omp_get_wtime() - dd;
		std::cout << "update time: " << dd << " sec" << std::endl;
		std::cout << "LSgrid matrices updated. t = " << t0 << " to " << (t0 + dt) << std::endl;
	}
}

//----------------------------------------------------------------------------------------------------------------------------------------


// SWD Routines


//----------------------------------------------------------------------------------------------------------------------------------------


#if defined(SWD_enable)

void sGrid::initialize_kinematics(SpectralWaveData *swd) {
	// Tell the swd object current application time...
	dx = (domain[1] - domain[0]) / std::max(1., double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));
	//dz = (domain_end[2] - domain_start[2]) / double(NZ - 1);
	ds = 1. / std::max(1., double(nl - 1));

	double dd = omp_get_wtime();

	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());

	// timestep T0.
	try {
		swd->UpdateTime(t0);
	}
	catch (SwdInputValueException& e) {  //Could be t > tmax from file.
		std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
		// If we will try again with a new value of t
		// we first need to call: swd.ExceptionClear()
		exit(EXIT_FAILURE);  // In this case we just abort.
	}

#pragma omp parallel // start parallel initialization
	{
		// timestep T0.
#pragma omp master
		std::cout << "Number of available threads: " << omp_get_num_threads() << std::endl;

		// Main grid
#pragma omp for collapse(2)
		for (int i = 0; i < nx; i++) {

			for (int j = 0; j < ny; j++) {
				double xpt = domain[0] + dx * i;
				double ypt = domain[2] + dy * j;
				double eta_temp = ETA0[i * ny + j] - swl;

				//std::cout << i << " " << j << ": " << ETA0[i * ny + j] << std::endl;


				for (int m = 0; m < nl; m++) {
					//std::cout << i << " " << j << " " << m << std::endl;
					double spt = s2tan(-1. + ds * m);
					double zpt = s2z(spt, eta_temp, water_depth);
					vector_swd U = swd->GradPhi(xpt, ypt, zpt);

					UX0[i * ny * nl + j * nl + m] = U.x;
					UY0[i * ny * nl + j * nl + m] = U.y;
					UZ0[i * ny * nl + j * nl + m] = U.z;

				}
			}
		}
	}

	// timestep T1.

	try {
		swd->UpdateTime(t0+dt);
	}
	catch (SwdInputValueException& e) {  //Could be t > tmax from file.
		std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
		// If we will try again with a new value of t
		// we first need to call: swd.ExceptionClear()
		exit(EXIT_FAILURE);  // In this case we just abort.
	}
#pragma omp parallel // start parallel initialization
	{
		// Main grid
#pragma omp for collapse(2)
		for (int i = 0; i < nx; i++) {

			for (int j = 0; j < ny; j++) {
				double xpt = domain[0] + dx * i;
				double ypt = domain[2] + dy * j;
				double eta_temp = ETA1[i * ny + j] - swl;


				for (int m = 0; m < nl; m++) {
					double spt = s2tan(-1. + ds * m);
					double zpt = s2z(spt, eta_temp, water_depth);
					vector_swd U = swd->GradPhi(xpt, ypt, zpt);

					UX1[i * ny * nl + j * nl + m] = U.x;
					UY1[i * ny * nl + j * nl + m] = U.y;
					UZ1[i * ny * nl + j * nl + m] = U.z;

				}
			}
		}
	} // End parallel initialization

	if (dump_vtk) {
		write_vtk(false);
		update_count++;
	}
	std::cout << "Generation of domain kinematics using SWD completed. ";
	dd = omp_get_wtime() - dd;
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
	std::cout << "Interpolation can commence..." << std::endl;
}

void sGrid::initialize_kinematics_with_ignore(SpectralWaveData* swd) {
	// Tell the swd object current application time...
	dx = (domain[1] - domain[0]) / std::max(1., double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));
	//dz = (domain_end[2] - domain_start[2]) / double(NZ - 1);
	ds = 1. / std::max(1., double(nl - 1));

	double dd = omp_get_wtime();

	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());

#pragma omp parallel // start parallel initialization
	{
		// timestep T0.
#pragma omp master
		std::cout << "Number of available threads: " << omp_get_num_threads() << std::endl;
		try {
			swd->UpdateTime(t0);
		}
		catch (SwdInputValueException& e) {  //Could be t > tmax from file.
			std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
			// If we will try again with a new value of t
			// we first need to call: swd.ExceptionClear()
			exit(EXIT_FAILURE);  // In this case we just abort.
		}

		// Main grid
#pragma omp for collapse(2)
		for (int i = 0; i < nx; i++) {

			for (int j = 0; j < ny; j++) {
				double xpt = domain[0] + dx * i;
				double ypt = domain[2] + dy * j;
				double eta_temp = ETA0[i * ny + j] - swl;

				if (!IGNORE[i * ny + j]) {
					for (int m = 0; m < nl; m++) {
						double spt = s2tan(-1. + ds * m);
						double zpt = s2z(spt, eta_temp, water_depth);
						vector_swd U = swd->GradPhi(xpt, ypt, zpt);

						UX0[i * ny * nl + j * nl + m] = U.x;
						UY0[i * ny * nl + j * nl + m] = U.y;
						UZ0[i * ny * nl + j * nl + m] = U.z;

					}
				}
				else {
					for (int m = 0; m < nl; m++) {
						UX0[i * ny * nl + j * nl + m] = 0.;
						UY0[i * ny * nl + j * nl + m] = 0.;
						UZ0[i * ny * nl + j * nl + m] = 0.;
					}
				}
			}
		}
		// timestep T1.
#pragma omp master
		try {
			swd->UpdateTime(t0 + dt);
		}
		catch (SwdInputValueException& e) {  //Could be t > tmax from file.
			std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
			// If we will try again with a new value of t
			// we first need to call: swd.ExceptionClear()
			exit(EXIT_FAILURE);  // In this case we just abort.
		}

		// Main grid
#pragma omp for collapse(2)
		for (int i = 0; i < nx; i++) {

			for (int j = 0; j < ny; j++) {
				double xpt = domain[0] + dx * i;
				double ypt = domain[2] + dy * j;
				double eta_temp = ETA1[i * ny + j] - swl;

				if (!IGNORE[i * ny + j]) {
					for (int m = 0; m < nl; m++) {
						double spt = s2tan(-1. + ds * m);
						double zpt = s2z(spt, eta_temp, water_depth);
						vector_swd U = swd->GradPhi(xpt, ypt, zpt);

						UX1[i * ny * nl + j * nl + m] = U.x;
						UY1[i * ny * nl + j * nl + m] = U.y;
						UZ1[i * ny * nl + j * nl + m] = U.z;

					}
				}
				else {
					for (int m = 0; m < nl; m++) {
						UX1[i * ny * nl + j * nl + m] = 0.;
						UY1[i * ny * nl + j * nl + m] = 0.;
						UZ1[i * ny * nl + j * nl + m] = 0.;
					}
				}
			}
		}
	} // End parallel initialization
	if (dump_vtk) {
		write_vtk(false);
		update_count++;
	}
	std::cout << "Generation of domain kinematics data completed. ";
	dd = omp_get_wtime() - dd;
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
	std::cout << "Interpolation can commence..." << std::endl;
}

void sGrid::initialize_surface_elevation(SpectralWaveData* swd, double t_target) {

	std::cout << "time: " << t_target << std::endl;
	t0 = t_target;

	// Allocating memory for storage of surface elevation and velocities

	dx = (domain[1] - domain[0]) / std::max(1., double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));

	bxmin = domain[0];
	bxmax = domain[1];
	bymin = domain[2];
	bymax = domain[3];

	double dd = omp_get_wtime();

	omp_set_num_threads(omp_get_max_threads());
	//omp_set_num_threads(1);

	// TIME T0
	try {
		swd->UpdateTime(t0);
	}
	catch (SwdInputValueException& e) {  //Could be t > tmax from file.
		std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
		// If we will try again with a new value of t
		// we first need to call: swd.ExceptionClear()
		exit(EXIT_FAILURE);  // In this case we just abort.
	}

#pragma omp parallel
	{
		// Main grid
#pragma omp for collapse(2) 
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				double xpt = domain[0] + dx * i;
				double ypt = domain[2] + dy * j;

				//std::cout << "wavelev: " << swd->Elev(xpt, ypt) << std::endl;

				ETA0[i * ny + j] = swd->Elev(xpt, ypt) + swl;
			}

		}
	}

	// Time T1
	try {
		swd->UpdateTime(t0 + dt);
	}
	catch (SwdInputValueException& e) {  //Could be t > tmax from file.
		std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
		// If we will try again with a new value of t
		// we first need to call: swd.ExceptionClear()
		exit(EXIT_FAILURE);  // In this case we just abort.
	}

#pragma omp parallel
	{
		// Main grid
#pragma omp for collapse(2)
		for (int i = 0; i < nx; i++) {

			for (int j = 0; j < ny; j++) {
				double xpt = domain[0] + dx * i;
				double ypt = domain[2] + dy * j;

				ETA1[i * ny + j] = swd->Elev(xpt, ypt) + swl;

			}
		}
	}
	dd = omp_get_wtime() - dd;

	std::cout << "Surface Elevation generated successfully. ";
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
}

void sGrid::initialize_surface_elevation_with_ignore(SpectralWaveData* swd, double t_target) {
	std::cout << "time: " << t_target << std::endl;
	t0 = t_target;

	// Allocating memory for storage of surface elevation and velocities

	dx = (domain[1] - domain[0]) / std::max(1., double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));

	bxmin = domain[0];
	bxmax = domain[1];
	bymin = domain[2];
	bymax = domain[3];

	double dd = omp_get_wtime();
	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());

#pragma omp parallel
	{
#pragma omp master
		try {
			swd->UpdateTime(t0);
		}
		catch (SwdInputValueException& e) {  //Could be t > tmax from file.
			std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
			// If we will try again with a new value of t
			// we first need to call: swd.ExceptionClear()
			exit(EXIT_FAILURE);  // In this case we just abort.
		}
		// Main grid
#pragma omp for collapse(2)
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (!IGNORE[i * ny + j]) {
					double xpt = domain[0] + dx * i;
					double ypt = domain[2] + dy * j;

					ETA0[i * ny + j] = swd->Elev(xpt, ypt) + swl;
				}
				else {
					ETA0[i * ny + j] = 0. + swl;
				}
			}
		}
#pragma omp master
		try {
			swd->UpdateTime(t0 + dt);
		}
		catch (SwdInputValueException& e) {  //Could be t > tmax from file.
			std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
			// If we will try again with a new value of t
			// we first need to call: swd.ExceptionClear()
			exit(EXIT_FAILURE);  // In this case we just abort.
		}
		// Main grid
#pragma omp for collapse(2)
		for (int i = 0; i < nx; i++) {

			for (int j = 0; j < ny; j++) {
				if (!IGNORE[i * ny + j]) {
					double xpt = domain[0] + dx * i;
					double ypt = domain[2] + dy * j;
					ETA1[i * ny + j] = swd->Elev(xpt, ypt) + swl;
				}
				else {
					ETA1[i * ny + j] = 0. + swl;
				}
			}
		}
	}
	dd = omp_get_wtime() - dd;

	std::cout << "Surface Elevation generated successfully. ";
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
}

void sGrid::update(SpectralWaveData* swd, double t_target)
{
	// Start by checking bounds
/*
if (!disable_checkbounds){
	CheckBounds();
}*/

// new time step
	if ((t_target / dt - (t0 + 2 * dt) / dt) > 0.) {
		double new_time = dt * std::floor(t_target / dt);
		initialize_surface_elevation(swd, new_time);
		initialize_kinematics(swd);
	}
	else {
		t0 += dt;
		// Updating surface elevations
		double dd = omp_get_wtime();
		//omp_set_num_threads(1);
		omp_set_num_threads(omp_get_max_threads());

		try {
			swd->UpdateTime(t0 + dt);
		}
		catch (SwdInputValueException& e) {  //Could be t > tmax from file.
			std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
			// If we will try again with a new value of t
			// we first need to call: swd.ExceptionClear()
			exit(EXIT_FAILURE);  // In this case we just abort.
		}
#pragma omp parallel
		{
			// Main grid
#pragma omp for collapse(2)
			for (int i = 0; i < nx; i++) {
				for (int j = 0; j < ny; j++) {
					if (!IGNORE[i * ny + j]) {
						double xpt = domain[0] + dx * i;
						double ypt = domain[2] + dy * j;
						ETA0[i * ny + j] = ETA1[i * ny + j];
						ETA1[i * ny + j] = swd->Elev(xpt, ypt) + swl;

						for (int m = 0; m < nl; m++) {
							double spt = s2tan(-1. + ds * m);
							double zpt = s2z(spt, ETA1[i * ny + j] - swl, water_depth);

							UX0[i * ny * nl + j * nl + m] = UX1[i * ny * nl + j * nl + m];
							UY0[i * ny * nl + j * nl + m] = UY1[i * ny * nl + j * nl + m];
							UZ0[i * ny * nl + j * nl + m] = UZ1[i * ny * nl + j * nl + m];

							vector_swd U = swd->GradPhi(xpt, ypt, zpt);

							UX1[i * ny * nl + j * nl + m] = U.x;
							UY1[i * ny * nl + j * nl + m] = U.y;
							UZ1[i * ny * nl + j * nl + m] = U.z;
						}
					}
				}
			}
		}

		if (dump_vtk) {
			write_vtk(false);
			update_count++;
		}

		dd = omp_get_wtime() - dd;
		std::cout << "update time: " << dd << " sec" << std::endl;
		std::cout << "LSgrid matrices updated. t = " << t0 << " to " << (t0 + dt) << std::endl;
	}
}

#endif



//----------------------------------------------------------------------------------------------------------------------------------------


// VTK output


//----------------------------------------------------------------------------------------------------------------------------------------


// when called, writes stored kinematics to file
void sGrid::write_vtk(bool endtime) {
	char buffer[256]; sprintf(buffer, "%05d", update_count);

	if (dirExists(vtk_directory_path.c_str()) == 0) {
		std::cout << "WARNING: Specified directory for storage of VTK files does not exist. Directory will be created at the following path:  " << vtk_directory_path << std::endl;
		createDirectory(vtk_directory_path);
	}


	std::string str(buffer);
	std::string fpath = (vtk_directory_path + vtk_prefix + buffer + ".vtu");
	std::cout << fpath << std::endl;
	FILE* fp = fopen(fpath.c_str(), "w");
	if (endtime)
		export_vtu(fp, true);
	else
		export_vtu(fp, false);
	fclose(fp);

	std::cout << "wrote kinematics to: " << fpath << std::endl;
}

/* exports sGrid at t= t0 to .vtu file for visualization in vtk/paraview */
void sGrid::export_vtu(FILE* fp, bool last)
{
	// write header
	fputs("<?xml version=\"1.0\"?>\n"
		"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
	fputs("\t <UnstructuredGrid>\n", fp);
	fprintf(fp, "\t\t <FieldData> \n");
	if (last) {
		fprintf(fp, "\t\t\t <DataArray type = \"Float64\" Name = \"%s\" NumberOfTuples = \"1\" format = \"ascii\" RangeMin = \"%.3f\" RangeMax = \"%.3f\"> \n", vtk_timelabel.c_str() , t0 + dt, t0 + dt);
		fprintf(fp, "\t\t\t %.3f \n", t0+dt);
	}
	else {
		fprintf(fp, "\t\t\t <DataArray type = \"Float64\" Name = \"%s\" NumberOfTuples = \"1\" format = \"ascii\" RangeMin = \"%.3f\" RangeMax = \"%.3f\"> \n", vtk_timelabel.c_str() , t0, t0);
		fprintf(fp, "\t\t\t %.3f \n", t0);
	}
	fprintf(fp, "\t\t\t </DataArray > \n");
	fprintf(fp, "\t\t </FieldData> \n");

	fprintf(fp, "\t\t <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nx*ny*nl, std::max((nx-1),1)* std::max((ny-1),1)* std::max((nl-1),1));
	
	// Loop over velocity data and store kinematics in cell vector stucture
	fputs("\t\t\t <PointData Scalars=\"scalars\">\n", fp);

	fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n");
	if (last) {
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				for (int m = 0; m < nl; m++) {
					fprintf(fp, "%g %g %g\n", UX1[i * ny * nl + j * nl + m], UY1[i * ny * nl + j * nl + m], UZ1[i * ny * nl + j * nl + m]);
				}
			}
		}
	}
	else {
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				for (int m = 0; m < nl; m++) {
					fprintf(fp, "%g %g %g\n", UX0[i * ny * nl + j * nl + m], UY0[i * ny * nl + j * nl + m], UZ0[i * ny * nl + j * nl + m]);
				}
			}
		}
	}
	fputs("\t\t\t\t </DataArray>\n", fp);

	fputs("\t\t\t </PointData>\n", fp);

	fputs("\t\t\t <Points>\n", fp);
	fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
	if (last) {
		for (int i = 0; i < nx; i++) {
			double xpt = domain[0] + dx * i;
			for (int j = 0; j < ny; j++) {
				double ypt = domain[2] + dy * j;
				double eta1_temp = ETA1[i * ny + j];
				for (int m = 0; m < nl; m++) {
					double spt = s2tan(-1. + ds * m);
					double zpt1 = s2z(spt, eta1_temp, water_depth);
					fprintf(fp, "%12.4f %12.4f %12.4f\n", xpt, ypt, zpt1);
				}
			}
		}
	}
	else {
		for (int i = 0; i < nx; i++) {
			double xpt = domain[0] + dx * i;
			for (int j = 0; j < ny; j++) {
				double ypt = domain[2] + dy * j;
				double eta0_temp = ETA0[i * ny + j];
				for (int m = 0; m < nl; m++) {
					double spt = s2tan(-1. + ds * m);
					double zpt0 = s2z(spt, eta0_temp, water_depth);
					fprintf(fp, "%12.4f %12.4f %12.4f\n", xpt, ypt, zpt0);
				}
			}
		}
	}
	fputs("\t\t\t\t </DataArray>\n", fp);
	fputs("\t\t\t </Points>\n", fp);

	fputs("\t\t\t <Cells>\n", fp);
	fputs("\t\t\t\t <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n", fp);

	if (nx > 1 && ny > 1 && nl > 1) {
		for (int i = 0; i < (nx - 1); i++) {
			for (int j = 0; j < (ny - 1); j++) {
				for (int m = 0; m < (nl - 1); m++) {
					int ape1 = nl * ny * i + nl * j + m;
					int ape2 = nl * ny * (i + 1) + nl * j + m;
					int ape3 = nl * ny * (i + 1) + nl * (j + 1) + m;
					int ape4 = nl * ny * i + nl * (j + 1) + m;
					int ape5 = nl * ny * i + nl * j + (m + 1);
					int ape6 = nl * ny * (i + 1) + nl * j + (m + 1);
					int ape7 = nl * ny * (i + 1) + nl * (j + 1) + (m + 1);
					int ape8 = nl * ny * i + nl * (j + 1) + (m + 1);
					fprintf(fp, "%u %u %u %u %u %u %u %u\n", ape1, ape2, ape3, ape4, ape5, ape6, ape7, ape8);
				}
			}
		}


		fputs("\t\t\t\t </DataArray>\n", fp);
		fputs("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);

		for (int i = 1; i < ((nx - 1) * (ny - 1) * (nl - 1) + 1); i++) {
			fprintf(fp, "%d \n", i * 8);
		}
		fputs("\t\t\t\t </DataArray>\n", fp);
		fputs("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
		for (int i = 1; i < ((nx - 1) * (ny - 1) * (nl - 1) + 1); i++) {
			fputs("12 \n", fp);
		}
	}
	// only single dimension i y direction.
	else if (nx > 1 && ny == 1 && nl > 1) {
		for (int i = 0; i < (nx - 1); i++) {		
			for (int m = 0; m < (nl - 1); m++) {
				int ape1 = nl * i + m;
				int ape2 = nl * (i + 1) + m;
				int ape3 = nl * (i + 1) + (m + 1);
				int ape4 = nl * i + (m + 1);
				fprintf(fp, "%u %u %u %u\n", ape1, ape2, ape3, ape4);
			}

		}

		fputs("\t\t\t\t </DataArray>\n", fp);
		fputs("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);

		for (int i = 1; i < ((nx - 1) * (nl - 1) + 1); i++) {
			fprintf(fp, "%d \n", i * 4);
		}
		fputs("\t\t\t\t </DataArray>\n", fp);
		fputs("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
		for (int i = 1; i < ((nx - 1) * (nl - 1) + 1); i++) {
			fputs("9 \n", fp);
		}

	}
	// only single dimension i x direction.
	else if (nx == 1 && ny > 1 && nl > 1) {
		for (int j = 0; j < (ny - 1); j++) {
			for (int m = 0; m < (nl - 1); m++) {
				int ape1 = nl * j + m;
				int ape2 = nl * (j + 1) + m;
				int ape3 = nl * (j + 1) + (m + 1);
				int ape4 = nl * j + (m + 1);
				fprintf(fp, "%u %u %u %u\n", ape1, ape2, ape3, ape4);
			}

		}

		fputs("\t\t\t\t </DataArray>\n", fp);
		fputs("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);

		for (int j = 1; j < ((ny - 1) * (nl - 1) + 1); j++) {
			fprintf(fp, "%d \n", j * 4);
		}
		fputs("\t\t\t\t </DataArray>\n", fp);
		fputs("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
		for (int j = 1; j < ((ny - 1) * (nl - 1) + 1); j++) {
			fputs("9 \n", fp);
		}
	}

	// only single dimension i z direction (lagrangian).
	else if (nx > 1 && ny > 1 && nl == 1) {
		for (int i = 0; i < (nx - 1); i++) {
			for (int j = 0; j < (ny - 1); j++) {
				int ape1 = ny * i + j;
				int ape2 = ny * (i + 1) + j;
				int ape3 = ny * (i + 1) + (j + 1);
				int ape4 = ny * i + (j + 1);
				fprintf(fp, "%u %u %u %u\n", ape1, ape2, ape3, ape4);
			}

		}

		fputs("\t\t\t\t </DataArray>\n", fp);
		fputs("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);

		for (int i = 1; i < ((nx - 1) * (ny - 1) + 1); i++) {
			fprintf(fp, "%d \n", i * 4);
		}
		fputs("\t\t\t\t </DataArray>\n", fp);
		fputs("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
		for (int i = 1; i < ((nx - 1) * (ny - 1) + 1); i++) {
			fputs("9 \n", fp);
		}

	}

	fputs("\t\t\t\t </DataArray>\n", fp);
	fputs("\t\t\t </Cells>\n", fp);
	fputs("\t\t </Piece>\n", fp);
	fputs("\t </UnstructuredGrid>\n", fp);
	fputs("</VTKFile>\n", fp);
	fflush(fp);
}

/* Set area of domain to ignore when update kinematics data. this is useful when prescribing kinematics at the boundaries*/
void sGrid::set_ignore()
{
	dx = (domain[1] - domain[0]) / std::max(1.,double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));

	for (int i = 0; i < nx; i++) {
		double xpt = domain[0] + dx * i;
		for (int j = 0; j < ny; j++) {
			double ypt = domain[2] + dy * j;
			//std::cout << xpt << " " << ypt << std::endl;
			//std::cout << domain_ignore[0] << " " << domain_ignore[1] << " " << domain_ignore[2] << " " << domain_ignore[3] << std::endl;
			if (xpt >= domain_ignore[0] && xpt <= domain_ignore[1] && ypt >= domain_ignore[2] && ypt <= domain_ignore[3]) {
				IGNORE[i * ny + j] = 1;
			}
			else {
				IGNORE[i * ny + j] = 0;
			}
		}
	}				
}

/******************************************************************************
 * Checks to see if a directory exists. Note: This method only checks the
 * existence of the full path AND if path leaf is a dir.
 *
 * @return  >0 if dir exists AND is a dir,
 *           0 if dir does not exist OR exists but not a dir,
 *          <0 if an error occurred (errno is also set)
 *****************************************************************************/
int sGrid::dirExists(const char* const path)
{
	struct stat info;

	int statRC = stat(path, &info);
	if (statRC != 0)
	{
		if (errno == ENOENT) { return 0; } // something along the path does not exist
		if (errno == ENOTDIR) { return 0; } // something in path prefix is not a dir
		return -1;
	}

	return (info.st_mode & S_IFDIR) ? 1 : 0;
}

void sGrid::createDirectory(std::string sPath) {

	
	int nError = 0;
#if defined(_WIN32)
	
	nError = _mkdir(sPath.c_str()); // can be used on Windows
#else
	mode_t nMode = 0733; // UNIX style permissions
	nError = mkdir(sPath.c_str(), nMode); // can be used on non-Windows
#endif
	if (nError != 0) {
		// handle your error here
	}
}