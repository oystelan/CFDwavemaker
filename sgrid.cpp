#include "sgrid.h"
#include "omp.h"
#include <algorithm>
#include <iostream>
#include <cmath>




// A streching function for setting variable layer thickness
double sGrid::slayer(int layerno) {
	//fprintf(stdout,"numlayers: %d",nl);
	double sfac = 3.0; // fixme: this should be made dimensionless and a function of specified wave
	double* dd = new double[NL];
	double ddsum = 0.;
	for (int ii = 0; ii < NL; ii++) {
		//fprintf(stdout,"%d",ii);
		dd[ii] = (1. + (NL - ii) * sfac);
		ddsum += dd[ii];
	}
	double layer_percentage = dd[layerno] / ddsum;
	delete[] dd;
	return layer_percentage;
}

// A function for constant layer thickness (equal spacing as a function of z)
double sGrid::clayer(int layerno) {
	//fprintf(stdout,"numlayers: %d",nl);
	return 1./double(NL);
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

// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 3
void sGrid::initialize_kinematics(Irregular& irregular) {
	// Allocating memory for storage of surface elevation and velocities
	UX0 = new double[NX * NY * NL];
	UY0 = new double[NX * NY * NL];
	UZ0 = new double[NX * NY * NL];

	UX1 = new double[NX * NY * NL];
	UY1 = new double[NX * NY * NL];
	UZ1 = new double[NX * NY * NL];

	

	std::cout << "Memory allocation successful for storage of kinematics." << std::endl;

	dx = (domain[1] - domain[0]) / double(NX - 1);
	dy = (domain[3] - domain[2]) / double(NY - 1);
	//dz = (domain_end[2] - domain_start[2]) / double(NZ - 1);
	ds = 1. / double(NL - 1);

	double dd = omp_get_wtime();

	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());

#pragma omp parallel // start parallel initialization
	{
#pragma omp master
		std::cout << "Number of available threads: " << omp_get_num_threads() << std::endl;

		// Main grid
#pragma omp for
		for (int i = 0; i < NX; i++) {
			double xpt = domain[0] + dx * i;
			for (int j = 0; j < NY; j++) {
				double ypt = domain[2] + dy * j;
				double eta0_temp = ETA0[i * NY + j];
				double eta1_temp = ETA1[i * NY + j];

				double Ux0 = irregular.u1(t0, xpt, ypt, 0.0) + irregular.u2(t0, xpt, ypt, 0.0);
				double Uy0 = irregular.v1(t0, xpt, ypt, 0.0) + irregular.v2(t0, xpt, ypt, 0.0);
				double Uz0 = irregular.w1(t0, xpt, ypt, 0.0) + irregular.w2(t0, xpt, ypt, 0.0);
				double Ux1 = irregular.u1(t0 + dt, xpt, ypt, 0.0) + irregular.u2(t0 + dt, xpt, ypt, 0.0);
				double Uy1 = irregular.v1(t0 + dt, xpt, ypt, 0.0) + irregular.v2(t0 + dt, xpt, ypt, 0.0);
				double Uz1 = irregular.w1(t0 + dt, xpt, ypt, 0.0) + irregular.w2(t0 + dt, xpt, ypt, 0.0);

				double PHI0_dxdz = irregular.phi1_dxdz(t0, xpt, ypt);
				double PHI0_dydz = irregular.phi1_dydz(t0, xpt, ypt);
				double PHI0_dzdz = irregular.phi1_dzdz(t0, xpt, ypt);

				double PHI1_dxdz = irregular.phi1_dxdz(t0 + dt, xpt, ypt);
				double PHI1_dydz = irregular.phi1_dydz(t0 + dt, xpt, ypt);
				double PHI1_dzdz = irregular.phi1_dzdz(t0 + dt, xpt, ypt);

				for (int m = 0; m < NL; m++) {
					double spt = s2tan(-1. + ds * m);
					double zpt0 = s2z(spt, eta0_temp, water_depth);
					double zpt1 = s2z(spt, eta1_temp, water_depth);

					if (zpt0 > 0.) {
						UX0[i * NY * NL + j * NL + m] = Ux0 + PHI0_dxdz * zpt0;
						UY0[i * NY * NL + j * NL + m] = Uy0 + PHI0_dydz * zpt0;
						UZ0[i * NY * NL + j * NL + m] = Uz0 + PHI0_dzdz * zpt0;
					}
					else {
						UX0[i * NY * NL + j * NL + m] = irregular.u1(t0, xpt, ypt, zpt0) + irregular.u2(t0, xpt, ypt, zpt0);
						UY0[i * NY * NL + j * NL + m] = irregular.v1(t0, xpt, ypt, zpt0) + irregular.v2(t0, xpt, ypt, zpt0);
						UZ0[i * NY * NL + j * NL + m] = irregular.w1(t0, xpt, ypt, zpt0) + irregular.w2(t0, xpt, ypt, zpt0);
					}
					if (zpt1 > 0.) {
						UX1[i * NY * NL + j * NL + m] = Ux1 + PHI1_dxdz * zpt1;
						UY1[i * NY * NL + j * NL + m] = Uy1 + PHI1_dydz * zpt1;
						UZ1[i * NY * NL + j * NL + m] = Uz1 + PHI1_dzdz * zpt1;
					}
					else {
						UX1[i * NY * NL + j * NL + m] = irregular.u1(t0 + dt, xpt, ypt, zpt1) + irregular.u2(t0 + dt, xpt, ypt, zpt1);
						UY1[i * NY * NL + j * NL + m] = irregular.v1(t0 + dt, xpt, ypt, zpt1) + irregular.v2(t0 + dt, xpt, ypt, zpt1);
						UZ1[i * NY * NL + j * NL + m] = irregular.w1(t0 + dt, xpt, ypt, zpt1) + irregular.w2(t0 + dt, xpt, ypt, zpt1);
					}
					
				}
			}
		}
	} // End parallel initialization

	std::cout << "Generation of domain kinematics data completed. ";
	dd = omp_get_wtime() - dd;
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
	std::cout << "Interpolation can commence..." << std::endl;
}


// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 3
void sGrid::initialize_surface_elevation(Irregular& irregular) {

	// Allocating memory for storage of surface elevation and velocities
	ETA0 = new double[NX * NY];
	ETA1 = new double[NX * NY];
	IGNORE = new int[NX * NY];
	std::cout << "Memory allocation successful for Surface elevation storage." << std::endl;

	dx = (domain[1] - domain[0]) / double(NX - 1);
	dy = (domain[3] - domain[2]) / double(NY - 1);

	double dd = omp_get_wtime();
	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());

#pragma omp parallel
	{
		// Main grid
#pragma omp for
		for (int i = 0; i < NX; i++) {
			double xpt = domain[0] + dx * i;
			for (int j = 0; j < NY; j++) {
				double ypt = domain[2] + dy * j;
				ETA0[i * NY + j] = irregular.eta1(t0, xpt, ypt) + irregular.eta2(t0, xpt, ypt);
				ETA1[i * NY + j] = irregular.eta1(t0+dt, xpt, ypt) + irregular.eta2(t0+dt, xpt, ypt);
				IGNORE[i * NY + j] = 0;
			}
		}
	}
	dd = omp_get_wtime() - dd;

	std::cout << "Surface Elevation generated successfully. ";
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
}

// When called, updates the arrays storing surface elevation and kinematics data for timestep t0 = t1, t1 = t1+dt
void sGrid::update(Irregular& irregular)
{
	// new time step
	t0 = t0 + dt;
	// Updating surface elevations
	double dd = omp_get_wtime();
	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel
	{
		// Main grid
#pragma omp for
		for (int i = 0; i < NX; i++) {
			double xpt = domain[0] + dx * i;
			for (int j = 0; j < NY; j++) {
				double ypt = domain[2] + dy * j;
				if (!IGNORE[i * NY + j]) {
					ETA0[i * NY + j] = ETA1[i * NY + j];
					ETA1[i * NY + j] = irregular.eta1(t0 + dt, xpt, ypt) + irregular.eta2(t0 + dt, xpt, ypt);

					double Ux1 = irregular.u1(t0 + dt, xpt, ypt, 0.0) + irregular.u2(t0 + dt, xpt, ypt, 0.0);
					double Uy1 = irregular.v1(t0 + dt, xpt, ypt, 0.0) + irregular.v2(t0 + dt, xpt, ypt, 0.0);
					double Uz1 = irregular.w1(t0 + dt, xpt, ypt, 0.0) + irregular.w2(t0 + dt, xpt, ypt, 0.0);

					double PHI1_dxdz = irregular.phi1_dxdz(t0 + dt, xpt, ypt);
					double PHI1_dydz = irregular.phi1_dydz(t0 + dt, xpt, ypt);
					double PHI1_dzdz = irregular.phi1_dzdz(t0 + dt, xpt, ypt);

					for (int m = 0; m < NL; m++) {
						double spt = s2tan(-1. + ds * m);
						double zpt1 = s2z(spt, ETA1[i * NY + j], water_depth);

						UX0[i * NY * NL + j * NL + m] = UX1[i * NY * NL + j * NL + m];
						UY0[i * NY * NL + j * NL + m] = UY1[i * NY * NL + j * NL + m];
						UZ0[i * NY * NL + j * NL + m] = UZ1[i * NY * NL + j * NL + m];

						if (zpt1 > 0.) {
							UX1[i * NY * NL + j * NL + m] = Ux1 + PHI1_dxdz * zpt1;
							UY1[i * NY * NL + j * NL + m] = Uy1 + PHI1_dydz * zpt1;
							UZ1[i * NY * NL + j * NL + m] = Uz1 + PHI1_dzdz * zpt1;
						}
						else {
							UX1[i * NY * NL + j * NL + m] = irregular.u1(t0 + dt, xpt, ypt, zpt1) + irregular.u2(t0 + dt, xpt, ypt, zpt1);
							UY1[i * NY * NL + j * NL + m] = irregular.v1(t0 + dt, xpt, ypt, zpt1) + irregular.v2(t0 + dt, xpt, ypt, zpt1);
							UZ1[i * NY * NL + j * NL + m] = irregular.w1(t0 + dt, xpt, ypt, zpt1) + irregular.w2(t0 + dt, xpt, ypt, zpt1);
						}
					}
				}
			}
		}
	}
	dd = omp_get_wtime() - dd;
	std::cout << "LSgrid matrices updated. t = " << t0 << " to " << (t0+dt) << std::endl;
}



/* Function for trilinear interpolation on a cartesian evenly spaced mesh*/
double sGrid::trilinear_interpolation(double* VAR, double xpt, double ypt, double zpt) {
	
	
	double wave_elev = eta(xpt, ypt);
	if (zpt > wave_elev) {
		return 0.;
	}
	double spt = tan2s(z2s(zpt, wave_elev, water_depth));
	double nxp = std::min(double(NX), std::max(0., (xpt - domain[0]) / dx));
	double nyp = std::min(double(NY), std::max(0., (ypt - domain[2]) / dy));
	double nsp =  (spt + 1.) / ds;

	double C000 = VAR[int(floor(nxp) * NY * NL + floor(nyp) * NL + floor(nsp))];
	double C001 = VAR[int(floor(nxp) * NY * NL + floor(nyp) * NL + ceil(nsp))];
	double C010 = VAR[int(floor(nxp) * NY * NL + ceil(nyp) * NL + floor(nsp))];
	double C011 = VAR[int(floor(nxp) * NY * NL + ceil(nyp) * NL + ceil(nsp))];
	double C100 = VAR[int(ceil(nxp) * NY * NL + floor(nyp) * NL + floor(nsp))];
	double C101 = VAR[int(ceil(nxp) * NY * NL + floor(nyp) * NL + ceil(nsp))];
	double C110 = VAR[int(ceil(nxp) * NY * NL + ceil(nyp) * NL + floor(nsp))];
	double C111 = VAR[int(ceil(nxp) * NY * NL + ceil(nyp) * NL + ceil(nsp))];
	double xd = nxp - floor(nxp);
	double yd = nyp - floor(nyp);
	double sd = nsp - floor(nsp);

	double C00 = C000 * (1. - xd) + C100 * xd;
	double C01 = C001 * (1. - xd) + C101 * xd;
	double C10 = C010 * (1. - xd) + C110 * xd;
	double C11 = C011 * (1. - xd) + C111 * xd;

	double C0 = C00 * (1. - yd) + C10 * yd;
	double C1 = C01 * (1. - yd) + C11 * yd;

	return C0 * (1. - sd) + C1 * sd;
}

/* bilinear interpolation function used to interpolate surface values on a regular evenly spaced grid*/
double sGrid::bilinear_interpolation(double* VAR, double xpt, double ypt) {

	double nxp = std::min(double(NX), std::max(0., (xpt - domain[0]) / dx));
	double nyp = std::min(double(NY), std::max(0., (ypt - domain[2]) / dy));

	double C00 = VAR[int(floor(nxp) * NY + floor(nyp))];
	double C01 = VAR[int(floor(nxp) * NY + ceil(nyp))];
	double C10 = VAR[int(ceil(nxp) * NY + floor(nyp))];
	double C11 = VAR[int(ceil(nxp) * NY + ceil(nyp))];

	double xd = nxp - floor(nxp);
	double yd = nyp - floor(nyp);

	double C0 = C00 * (1. - xd) + C10 * xd;
	double C1 = C01 * (1. - xd) + C11 * xd;

	return C0 * (1. - yd) + C1 * yd;
}


double sGrid::u(double xpt, double ypt, double zpt) {
	return trilinear_interpolation(UX0, xpt, ypt, zpt);
}

double sGrid::v(double xpt, double ypt, double zpt) {
	return trilinear_interpolation(UY0, xpt, ypt, zpt);
}

double sGrid::w(double xpt, double ypt, double zpt) {
	return trilinear_interpolation(UZ0, xpt, ypt, zpt);
}

double sGrid::eta(double xpt, double ypt) {
	return bilinear_interpolation(ETA0, xpt, ypt);
}


bool sGrid::CheckTime(double tpt) {
	/* Checks to see if the time tpt is within the interval t0 to t1. If so, returns true*/
	if (tpt > t0+dt) {
		std::cout << "T0: " << t0 << ", t1: " << (t0 + dt) << ", tpt: " << tpt << std::endl;
		return false;
		
	}
	return true;
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

/* Set area of domain to ignore when update kinematics data. this is useful when prescribing kinematics at the boundaries*/
void sGrid::set_ignore(double* bounds)
{
	for (int i = 0; i < NX; i++) {
		double xpt = domain[0] + dx * i;
		for (int j = 0; j < NY; j++) {
			double ypt = domain[2] + dy * j;
			if (xpt >= bounds[0] && xpt <= bounds[1] && ypt >= bounds[2] && ypt <= bounds[3]) {
				IGNORE[i * NY + j] = 1;
			}
		}
	}				
}