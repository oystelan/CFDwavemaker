
#include "Utils.h"



#include "omp.h"
#include <algorithm>
#include <iostream>
#include <cmath>

/***
void Grid::update_boundary_arrays(Irregular &irregular, double tpt) {
	t0 = t1;
	t1 += dt;

}

// If point which is outside the domain of the bounding box is asked for, the bounding box is redifined and reinitialized
void Grid::redefine_boundary_wallx(Irregular &irregular, double tpt, double xpt, double ypt, double zpt)
{
	// adjust wallx array size to include the new point
	wallxsize[0] = std::min(wallxsize[0], xpt);
	wallxsize[1] = std::max(wallxsize[1], xpt);
	wallxsize[2] = std::min(wallxsize[2], ypt);
	wallxsize[3] = std::max(wallxsize[3], ypt);
	wallxsize[4] = std::min(wallxsize[4], zpt);
	wallxsize[5] = std::max(wallxsize[5], zpt);

	
	init_boundary_wallx(irregular, tpt);
}

void Grid::update_boundary_wallx(Irregular &irregular, double tpt) {
	// updating timesteps
	t0 = t1;
	t1 = t1 + dt;

	double dd = omp_get_wtime();

	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());

	#pragma omp parallel // start parallel initialization
	{
		//#pragma omp master
		//std::cout << "Number of available threads: " << omp_get_num_threads() << std::endl;

		double xpt, ypt, zpt;
		double eta0_temp, eta1_temp;

		// Main grid
		#pragma omp for
		for (int i = 0; i < wallx_nx; i++) {
			xpt = wallxsize[0] + dx_wx * i;
			for (int j = 0; j < wallx_ny; j++) {
				ypt = wallxsize[2] + dy_wx * j;

				eta1_temp = irregular.eta(t1, xpt, ypt);

				double Ux1 = irregular.u1(t1, xpt, ypt, 0.0) + irregular.u2(t1, xpt, ypt, 0.0);
				double Uy1 = irregular.v1(t1, xpt, ypt, 0.0) + irregular.v2(t1, xpt, ypt, 0.0);
				double Uz1 = irregular.w1(t1, xpt, ypt, 0.0) + irregular.w2(t1, xpt, ypt, 0.0);

				double PHI1_dxdz = irregular.phi1_dxdz(t1, xpt, ypt);
				double PHI1_dydz = irregular.phi1_dydz(t1, xpt, ypt);
				double PHI1_dzdz = irregular.phi1_dzdz(t1, xpt, ypt);

				ETAX0[i * wallx_ny + j] = ETAX1[i * wallx_ny + j];

				ETAX1[i * wallx_ny + j] = eta1_temp;

				for (int m = 0; m < NZ; m++) {
					UX0[i * wallx_ny * wallx_nz + j * wallx_nz + m] = UX1[i * wallx_ny * wallx_nz + j * wallx_nz + m];
					VX0[i * wallx_ny * wallx_nz + j * wallx_nz + m] = VX1[i * wallx_ny * wallx_nz + j * wallx_nz + m];
					WX0[i * wallx_ny * wallx_nz + j * wallx_nz + m] = WX1[i * wallx_ny * wallx_nz + j * wallx_nz + m];

					zpt = wallxsize[4] + dz_wx * m;
					if (zpt > (std::max(ETAX0[i * wallx_ny + j], eta1_temp) + dz_wx)) {
						UX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = 0.0;
						VX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = 0.0;
						WX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = 0.0;
					}
					else if (zpt > 0.) {
						UX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = Ux1 + PHI1_dxdz * zpt;
						VX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = Uy1 + PHI1_dydz * zpt;
						WX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = Uz1 + PHI1_dzdz * zpt;
					}
					else {
						UX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = irregular.u1(t1, xpt, ypt, zpt) + irregular.u2(t1, xpt, ypt, zpt);
						VX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = irregular.v1(t1, xpt, ypt, zpt) + irregular.v2(t1, xpt, ypt, zpt);
						WX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = irregular.w1(t1, xpt, ypt, zpt) + irregular.w2(t1, xpt, ypt, zpt);
					}
				}
			}
		}
	} // End parallel initialization
	dd = omp_get_wtime() - dd;
	std::cout << "updated wallX boundaries in " << dd << " seconds." << std::endl;
}


void Grid::allocate_wallx_memory() {
	// arrays for time t= t0
	UX0 = new double[wallx_nx * wallx_ny * wallx_nz];
	VX0 = new double[wallx_nx * wallx_ny * wallx_nz];
	WX0 = new double[wallx_nx * wallx_ny * wallx_nz];
	ETAX0 = new double [wallx_nx * wallx_ny];
	for (int i = 0; i < (wallx_nx * wallx_ny); i++) {
		ETAX0[i] = double(i);
	}
	// arrays for time t= t1
	UX1 = new double[wallx_nx * wallx_ny * wallx_nz];
	VX1 = new double[wallx_nx * wallx_ny * wallx_nz];
	WX1 = new double[wallx_nx * wallx_ny * wallx_nz];
	ETAX1 = new double[wallx_nx * wallx_ny];

	std::cout << "WallX memory allocated" << std::endl;
}

void Grid::init_boundary_wallx(Irregular &irregular, double tpt) {
	t0 = tpt;
	t1 = tpt + dt;

	dx_wx = (wallxsize[1] - wallxsize[0]) / double(wallx_nx - 1);
	dy_wx = (wallxsize[3] - wallxsize[2]) / double(wallx_ny - 1);
	dz_wx = (wallxsize[5] - wallxsize[4]) / double(wallx_nz - 1);
	
	double dd = omp_get_wtime();

	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());

	#pragma omp parallel // start parallel initialization
	{
		#pragma omp master
		std::cout << "Number of available threads: " << omp_get_num_threads() << std::endl;

		double xpt, ypt, zpt;
		double eta0_temp, eta1_temp;

		// Main grid
		#pragma omp for
		for (int i = 0; i < wallx_nx; i++) {
			xpt = wallxsize[0] + dx_wx * i;
			for (int j = 0; j < wallx_ny; j++) {
				ypt = wallxsize[2] + dy_wx * j;
				eta0_temp = irregular.eta(t0, xpt, ypt);
				eta1_temp = irregular.eta(t1, xpt, ypt);			
				
				double Ux0 = irregular.u1(t0, xpt, ypt, 0.0) + irregular.u2(t0, xpt, ypt, 0.0);
				double Uy0 = irregular.v1(t0, xpt, ypt, 0.0) + irregular.v2(t0, xpt, ypt, 0.0);
				double Uz0 = irregular.w1(t0, xpt, ypt, 0.0) + irregular.w2(t0, xpt, ypt, 0.0);
				double Ux1 = irregular.u1(t1, xpt, ypt, 0.0) + irregular.u2(t1, xpt, ypt, 0.0);
				double Uy1 = irregular.v1(t1, xpt, ypt, 0.0) + irregular.v2(t1, xpt, ypt, 0.0);
				double Uz1 = irregular.w1(t1, xpt, ypt, 0.0) + irregular.w2(t1, xpt, ypt, 0.0);

				double PHI0_dxdz = irregular.phi1_dxdz(t0, xpt, ypt);
				double PHI0_dydz = irregular.phi1_dydz(t0, xpt, ypt);
				double PHI0_dzdz = irregular.phi1_dzdz(t0, xpt, ypt);
				double PHI1_dxdz = irregular.phi1_dxdz(t1, xpt, ypt);
				double PHI1_dydz = irregular.phi1_dydz(t1, xpt, ypt);
				double PHI1_dzdz = irregular.phi1_dzdz(t1, xpt, ypt);
				
				//std::cout << ETAX0[i * wallx_ny + j] << std::endl;
				ETAX0[i * wallx_ny + j] = eta0_temp;
				ETAX1[i * wallx_ny + j] = eta1_temp;

				
				
				for (int m = 0; m < NZ; m++) {
					zpt = wallxsize[4] + dz_wx * m;
					if (zpt > (std::max(eta0_temp,eta1_temp) + dz_wx)) {
						UX0[i * wallx_ny * wallx_nz + j * wallx_nz + m] = 0.0;
						VX0[i * wallx_ny * wallx_nz + j * wallx_nz + m] = 0.0;
						WX0[i * wallx_ny * wallx_nz + j * wallx_nz + m] = 0.0;
						UX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = 0.0;
						VX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = 0.0;
						WX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = 0.0;
					}
					else if (zpt > 0.) {
						UX0[i * wallx_ny * wallx_nz + j * wallx_nz + m] = Ux0 + PHI0_dxdz * zpt;
						VX0[i * wallx_ny * wallx_nz + j * wallx_nz + m] = Uy0 + PHI0_dydz * zpt;
						WX0[i * wallx_ny * wallx_nz + j * wallx_nz + m] = Uz0 + PHI0_dzdz * zpt;
						UX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = Ux1 + PHI1_dxdz * zpt;
						VX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = Uy1 + PHI1_dydz * zpt;
						WX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = Uz1 + PHI1_dzdz * zpt;
					}
					else {
						UX0[i * wallx_ny * wallx_nz + j * wallx_nz + m] = irregular.u1(t0, xpt, ypt, zpt) + irregular.u2(t0, xpt, ypt, zpt);
						VX0[i * wallx_ny * wallx_nz + j * wallx_nz + m] = irregular.v1(t0, xpt, ypt, zpt) + irregular.v2(t0, xpt, ypt, zpt);
						WX0[i * wallx_ny * wallx_nz + j * wallx_nz + m] = irregular.w1(t0, xpt, ypt, zpt) + irregular.w2(t0, xpt, ypt, zpt);
						UX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = irregular.u1(t1, xpt, ypt, zpt) + irregular.u2(t1, xpt, ypt, zpt);
						VX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = irregular.v1(t1, xpt, ypt, zpt) + irregular.v2(t1, xpt, ypt, zpt);
						WX1[i * wallx_ny * wallx_nz + j * wallx_nz + m] = irregular.w1(t1, xpt, ypt, zpt) + irregular.w2(t1, xpt, ypt, zpt);
					}
				}
			}
		}
	} // End parallel initialization
	dd = omp_get_wtime() - dd;
	std::cout << "Initialized wallX boundaries in " << dd << " seconds." << std::endl;
}



// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 3
void Grid::initialize_kinematics(Irregular &irregular, double tpt) {
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
						UX[i * NY * NZ + j * NZ + m] = irregular.u1(tpt, xpt, ypt, zpt) + irregular.u2(tpt, xpt, ypt, zpt);
						UY[i * NY * NZ + j * NZ + m] = irregular.v1(tpt, xpt, ypt, zpt) + irregular.v2(tpt, xpt, ypt, zpt);
						UZ[i * NY * NZ + j * NZ + m] = irregular.w1(tpt, xpt, ypt, zpt) + irregular.w2(tpt, xpt, ypt, zpt);
					}
					//UX[i*NY*NZ + j*NZ + m] = uu(tpt, xpt, ypt, zpt);
					//UY[i*NY*NZ + j*NZ + m] = vv(tpt, xpt, ypt, zpt);
					//UZ[i*NY*NZ + j*NZ + m] = ww(tpt, xpt, ypt, zpt);
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
					UXL[i * NYL * NZL + j * NZL + m] = irregular.u1(tpt, xpt, ypt, zpt) + irregular.u2(tpt, xpt, ypt, zpt);
					UYL[i * NYL * NZL + j * NZL + m] = irregular.v1(tpt, xpt, ypt, zpt) + irregular.v2(tpt, xpt, ypt, zpt);
					UZL[i * NYL * NZL + j * NZL + m] = irregular.w1(tpt, xpt, ypt, zpt) + irregular.w2(tpt, xpt, ypt, zpt);
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
void Grid::initialize_surface_elevation(Irregular &irregular, double tpt) {

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
				ETA[i * NY + j] = irregular.eta1(tpt, xpt, ypt) + irregular.eta2(tpt, xpt, ypt);
			}
		}
	}
	dd = omp_get_wtime() - dd;

	std::cout << "Surface Elevation generated successfully. ";
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
	initsurf = 1;
}

// Function for trilinear interpolation on a cartesian evenly spaced mesh
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


// Function for trilinear interpolation on a cartesian evenly spaced mesh on the lower part of the domain
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
}



// bilinear interpolation function used to interpolate surface values on a regular evenly spaced grid
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
	return bilinear_interpolation(ETA, xpt, ypt, NX, NY, dx, dy, domain_start);
}

double Grid::eta_wall(double tpt, double xpt, double ypt) {
	double eta0 = bilinear_interpolation(ETAX0, xpt, ypt, wallx_nx, wallx_ny, dx_wx, dy_wx, wallx_start);
	double eta1 = bilinear_interpolation(ETAX1, xpt, ypt, wallx_nx, wallx_ny, dx_wx, dy_wx, wallx_start);

	double yd = (tpt - t0) / (t1 - t0);
	//std::cout << eta0 << " " << eta1 << " " << yd << std::endl;
	return  eta0 * (1. - yd) + eta1 * yd;
}

double Grid::u_wall(double tpt, double xpt, double ypt, double zpt) {
	double u0 = trilinear_interpolation(UX0, xpt, ypt, zpt, wallx_nx, wallx_ny, wallx_nz, dx_wx, dy_wx, dz_wx, wallx_start);
	double u1 = trilinear_interpolation(UX1, xpt, ypt, zpt, wallx_nx, wallx_ny, wallx_nz, dx_wx, dy_wx, dz_wx, wallx_start);
	double yd = (tpt - t0) / (t1 - t0);
	return  u0 * (1. - yd) + u1 * yd;
}

double Grid::v_wall(double tpt, double xpt, double ypt, double zpt) {
	double v0 = trilinear_interpolation(VX0, xpt, ypt, zpt, wallx_nx, wallx_ny, wallx_nz, dx_wx, dy_wx, dz_wx, wallx_start);
	double v1 = trilinear_interpolation(VX1, xpt, ypt, zpt, wallx_nx, wallx_ny, wallx_nz, dx_wx, dy_wx, dz_wx, wallx_start);
	double yd = (tpt - t0) / (t1 - t0);
	return  v0 * (1. - yd) + v1 * yd;
}

double Grid::w_wall(double tpt, double xpt, double ypt, double zpt) {
	double w0 = trilinear_interpolation(WX0, xpt, ypt, zpt, wallx_nx, wallx_ny, wallx_nz, dx_wx, dy_wx, dz_wx, wallx_start);
	double w1 = trilinear_interpolation(WX1, xpt, ypt, zpt, wallx_nx, wallx_ny, wallx_nz, dx_wx, dy_wx, dz_wx, wallx_start);
	double yd = (tpt - t0) / (t1 - t0);
	return  w0 * (1. - yd) + w1 * yd;
}


bool Grid::CheckTime(double tpt) {
	// Checks to see if the time tpt is within the interval t0 to t1. If so, returns true
	if (tpt > t1) {
		return false;
		std::cout << "T0: " << t0 << ", t1: " << t1 << ", tpt: " << tpt << std::endl;
	}
	else{
		return true;
	}
}

// function to find if given point 
// lies inside a given rectangle or not. 
bool Grid::CheckBounds(double* bounds, double x, double y, double z)
{
	if (x >= bounds[0] && x <= bounds[1] && y >= bounds[2] && y <= bounds[3] && z >= bounds[4] && z <= bounds[5])
		return true;
	else {
		std::cout << "position x,y,z=" << x << "," << y << "," << z << " out of bounds. updating wallx borders" << std::endl;
		return false;
	}
}

*/

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
		return *std::min_element(std::begin(ramps), std::end(ramps));
	}
	else {
		return 1.0;
	}
}
