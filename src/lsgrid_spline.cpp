#include "lsgrid_spline.h"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <string>
#include "Utils.h"
#include "omp.h"


//----------------------------------------------------------------------------------------------------------------------------------------


// LSgrid functions


//----------------------------------------------------------------------------------------------------------------------------------------



// Allocation of memory to storage matrices
void lsGrid::allocate() {
	// Primary fields
	ETA = new double[4 * nx * ny];
	UX = new double[4 * nx * ny * nl];
	UY = new double[4 * nx * ny * nl];
	UZ = new double[4 * nx * ny * nl];

	// secondary fields (derivatives)
	etadx = new double[2 * nx * ny];
	etady = new double[2 * nx * ny];
	etadt = new double[4 * nx * ny];
	udt = new double[4 * nx * ny * nl];
	udx = new double[2 * nx * ny * nl];
	udy = new double[2 * nx * ny * nl];
	uds = new double[2 * nx * ny * nl];
	vdt = new double[4 * nx * ny * nl];
	vdx = new double[2 * nx * ny * nl];
	vdy = new double[2 * nx * ny * nl];
	vds = new double[2 * nx * ny * nl];
	wdt = new double[4 * nx * ny * nl];
	wdx = new double[2 * nx * ny * nl];
	wdy = new double[2 * nx * ny * nl];
	wds = new double[2 * nx * ny * nl];

	// Misc fields
	IGNORE = new int[nx * ny];
}

// A streching function for setting variable layer thickness
double lsGrid::slayer(int layerno) {
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
double lsGrid::clayer(int layerno) {
	//fprintf(stdout,"numlayers: %d",nl);
	return 1./double(nl);
}

// Transforms normal z axis to streched sigma coordinates 
// s defined between -1 (seabed) and 0 (free surface)
double lsGrid::z2s(double z, double wave_elev, double depth) {
	return (z - (wave_elev - swl)) / (depth + (wave_elev - swl));
}

// Transforms stretched sigma coordinate to normal z
double lsGrid::s2z(double s, double wave_elev, double depth) {
	return (wave_elev - swl) + s * (depth + (wave_elev - swl));
}

double lsGrid::s2tan(double s) {
	// s defined between -1 and 0, where -1 is seabed, 0 is sea surface
	// returns tangens strethced coordintates tan, which is also defined between -1 and 0	
	return -std::pow(std::tan(-s * tan_a) , tan_b) / std::pow(std::tan(tan_a), tan_b);
}

double lsGrid::tan2s(double t) {
	// The inverse of the above function s2tan. from tan stretched to normal constant spacing	
	return -std::atan(std::pow(-t * std::pow(std::tan(tan_a) , tan_b) , 1. / tan_b)) / tan_a;
}

double lsGrid::cart2sigmaS(double zpt, double wave_elev, double depth){
    double z = std::min(zpt, wave_elev); // impose an upper bound
    double sigma = std::max(z2s(z, wave_elev, depth), -1); // transform to sigma coordinates
    return 1 - std::atan(std::pow(-sigma * std::pow(std::tan(tan_a), tan_b), 1. / tan_b)) / tan_a;
    // defined now from 0 to 1 where 0 is sea bed and 1 is surface
}

double lsGrid::sigmaS2cart(double s, double wave_elev, double depth){
    double sigma = -std::pow(std::tan(-(s - 1) * tan_a) , tan_b) / std::pow(std::tan(tan_a), tan_b);
	return s2z(sigma, wave_elev, depth);
}

void lsGrid::update_gradient_dt(double* data, double* graddt){
    // function for computation of temporal gradients using central difference
    if (update_count > 0){
        int tid = (update_count + 2) % 4;
        int tid1 = (update_count + 3) % 4;
        int tid0 = (update_count + 1) % 4;
		for (int i = 0; i < nx*ny*nz; i++) {
            graddt[tid * nx * ny * nz + i] = (data[tid1 * nx * ny * nz + i] - data[tid0 * nx * ny * nz + i]) / (2 * dt);
		}
	}
    else {
        for (int i = 0; i < nx*ny*nz; i++) {
            graddt[1 * nx * ny * nz + i] = (data[2 * nx * ny * nz + i] - data[1 * nx * ny * nz + i]) / dt;
            graddt[2 * nx * ny * nz + i] = (data[3 * nx * ny * nz + i] - data[1 * nx * ny * nz + i]) / (2 * dt);
		}
	}
}

void lsGrid::update_gradient_eta_dt(double* data, double* graddt){
    // function for computation of temporal gradients using central difference
    if (update_count > 0){
        int tid = (update_count + 2) % 4;
        int tid1 = (update_count + 3) % 4;
        int tid0 = (update_count + 1) % 4;
        for (int i = 0; i < nx*ny; i++) {
            graddt[tid * nx * ny + i] = (data[tid1 * nx * ny + i] - data[tid0 * nx * ny + i]) / (2 * dt);
		}
	}
    else {
        for (int i = 0; i < nx*ny; i++) {
            graddt[1 * nx * ny + i] = (data[2 * nx * ny + i] - data[1 * nx * ny + i]) / dt;
            graddt[2 * nx * ny + i] = (data[3 * nx * ny + i] - data[1 * nx * ny + i]) / (2 * dt);
		}
	}
}

void lsGrid::update_gradient_eta_dxdy(double* eta, double* gradx,double* grady){
    // tid: specifies which time step to compute gradients for (value from 0 to 1)
    // compute spatial gradient of the sea surface using central difference

    int tid0 = (update_count + 1) % 4;
    int tid1 = (update_count + 2) % 4;
	for (int i = 1; i < (nx-1); i++) {
		for (int j = 0; j < ny; j++) {
            gradx[0 * nx * ny + i * ny + j] = (eta[tid0 * nx * ny + (i + 1) * ny + j] -
                                                   eta[tid0 * nx * ny +(i - 1) * ny + j]) / (2 * dx);
            gradx[1 * nx * ny + i * ny + j] = (eta[tid1 * nx * ny + (i + 1) * ny + j] -
                                                 eta[tid1 * nx * ny + (i - 1) * ny + j]) / (2 * dx);
		}
	}
    
	// boundaries are special cases (x=0 and x=xn)
    for (int j = 0; j < ny; j++) {
        gradx[0 * nx * ny + j] = (eta[tid0 * nx * ny + ny + j] - eta[tid0 * nx * ny +j]) / dx;
        gradx[0 * nx * ny + (nx-1) * ny + j] = (eta[tid0 * nx * ny +(nx-1) * ny + j] -
                                                   eta[tid0 * nx * ny +(nx-2) * ny + j]) / dx;
        gradx[1 * nx * ny + j] = (eta[tid1 * nx * ny + ny + j] - eta[tid1 * nx * ny + j]) / dx;
        gradx[1 * nx * ny + (nx - 1) * ny + j] = (eta[tid1 * nx * ny + (nx - 1) * ny + j] -
                                                   eta[tid1 * nx * ny + (nx - 2) * ny + j]) / dx;
	}

    // compute spatial gradient using central difference
    for (int i = 0; i < nx; i++) {
		for (int j = 1; j < (ny-1); j++) {
            grady[0 * nx * ny +i * ny + j] = (eta[tid0 * nx * ny + i * ny + (j + 1)] -
                                                  eta[tid0 * nx * ny + i * ny + (j - 1)]) / (2 * dy)
            grady[1 * nx * ny + i * ny + j] = (eta[tid1 * nx * ny + i * ny + (j + 1)] -
                                                 eta[tid1 * nx * ny + i * ny + (j - 1)]) / (2 * dy)
		}
	}
    // y=0 and y=yn
    for (int i = 0; i < nx; i++) {
        grady[0 * nx * ny + i * ny] = (eta[tid0 * nx * ny + i * ny + 1] - eta[tid0 * nx * ny + i * ny]) / dy
        grady[0 * nx * ny + i * ny + (ny - 1)] = (eta[tid0 * nx * ny + i * ny + (ny - 1)] -
                                                     eta[tid0 * nx * ny +i * ny + (ny - 2)]) / dy
        grady[1 * nx * ny + i * ny] = (eta[tid1 * nx * ny + i * ny + 1] - eta[tid1 * nx * ny + i * ny]) / dy
        grady[1 * nx * ny + i * ny + (ny - 1)] = (eta[tid1 * nx * ny + i * ny + (ny - 1)] -
                                                   eta[tid1 * nx * ny + i * ny + (ny - 2)]) / dy

	}
}

void lsGrid::update_gradient_dxdydz(double* data, double* gradx, double* grady, double* gradz){
    // compute spatial gradient using central difference

    int tid0 = (update_count + 1) % 4;
    int tid1 = (update_count + 2) % 4;
    for (int i = 1; i < (nx-1); i++) {
		for (int j = 0; j < ny; j++) {
			for (int m = 0; m < nl; m++) {
                gradx[0 * nx * ny * nl + i * ny * nl + j * nl + m] = ((data[tid0 * nx * ny * nl + (i + 1) * ny * nl + j * nl + m] -
                    data[tid0 * nx * ny * nl + (i - 1) * ny * nl + j * nl + m]) / (2 * dx));
                gradx[1 * nx * ny * nl + i * ny * nl + j * nl + m] = ((data[tid1 * nx * ny * nl + (i + 1) * ny * nl + j * nl + m] -
                      data[tid1 * nx * ny * nl + (i - 1) * ny * nl + j * nl + m]) / (2 * dx));
			}
		}
	}
    // x=0 and x=xn
    for (int j = 0; j < ny; j++) {
		for (int m = 0; m < nl; m++) {
            gradx[0 * nx * ny * nl + j * nl + m] = (
                data[tid0 * nx * ny * nl + ny * nl + j * nl + m] -
                data[tid0 * nx * ny * nl + j * nl + m]) / dx;
            gradx[0 * nx * ny * nl + (nx-1) * ny * nl + j * nl + m] = (
                data[tid0 * nx * ny * nl + (nx-1) * ny * nl + j * nl + m] -
                data[tid0 * nx * ny * nl + (nx-2) * ny * nl + j * nl + m]) / dx;
            gradx[1 * nx * ny * nl + j * nl + m] = (
                 data[tid1 * nx * ny * nl + ny * nl + j * nl + m] -
                 data[tid1 * nx * ny * nl + j * nl + m]) / dx;
            gradx[1 * nx * ny * nl + (nx - 1) * ny * nl + j * nl + m] = (
                  data[tid1 * nx * ny * nl + (nx - 1) * ny * nl + j * nl + m] -
                  data[tid1 * nx * ny * nl + (nx - 2) * ny * nl + j * nl + m]) / dx;
		}
	}

    // compute spatial gradient using central difference
    for (int i = 0; i < nx; i++) {
		for (int j = 1; j < (ny-1); j++) {
			for (int m = 0; m < nl; m++) {
                grady[0 * nx * ny * nl + i * ny * nl + j * nl + m] = (
                        (data[tid0 * nx * ny * nl + i * ny * nl + (j + 1) * nl + m] -
                         data[tid0 * nx * ny * nl + i * ny * nl + (j - 1) * nl + m]) / (2 * dy));
                grady[1 * nx * ny * nl + i * ny * nl + j * nl + m] = (
                        (data[tid1 * nx * ny * nl + i * ny * nl + (j + 1) * nl + m] -
                         data[tid1 * nx * ny * nl + i * ny * nl + (j - 1) * nl + m]) / (2 * dy));
			}
		}
	}

    // y=0 and y=yn
    for (int i = 0; i < nx; i++) {
		for (int m = 0; m < nl; m++) {
            grady[0 * nx * ny * nl + i * ny * nl + m] = (
                    (data[tid0 * nx * ny * nl + i * ny * nl + nl + m] -
                     data[tid0 * nx * ny * nl + i * ny * nl + m]) / dy);
            grady[0 * nx * ny * nl + i * ny * nl + (ny - 1) * nl + m] = (
                    (data[tid0 * nx * ny * nl + i * ny * nl + (ny - 1) * nl + m] -
                     data[tid0 * nx * ny * nl + i * ny * nl + (ny - 2) * nl + m]) / dy);
            grady[1 * nx * ny * nl + i * ny * nl + m] = (
                    (data[tid1 * nx * ny * nl + i * ny * nl + nl + m] -
                     data[tid1 * nx * ny * nl + i * ny * nl + m]) / dy);
            grady[1 * nx * ny * nl + i * ny * nl + (ny - 1) * nl + m] = (
                    (data[tid1 * nx * ny * nl + i * ny * nl + (ny - 1) * nl + m] -
                     data[tid1 * nx * ny * nl + i * ny * nl + (ny - 2) * nl + m]) / dy);
		}
	}
    
	// compute grad z spatial gradient using central difference
     for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int m = 1; m < (nl-1); m++) {
                gradz[0 * nx * ny * nl + i * ny * nl + j * nl + m] = (
                        (data[tid0 * nx * ny * nl + i * ny * nl + j * nl + (m + 1)] -
                         data[tid0 * nx * ny * nl + i * ny * nl + j * nl + (m - 1)]) / (2 * ds));
                gradz[1 * nx * ny * nl + i * ny * nl + j * nl + m] = (
                        (data[tid1 * nx * ny * nl + i * ny * nl + j * nl + (m + 1)] -
                         data[tid1 * nx * ny * nl + i * ny * nl + j * nl + (m - 1)]) / (2 * ds));
			}
		}
	 }
    // z=0 and z=zn
    for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
            gradz[0 * nx * ny * nl + i * ny * nl + j * nl] = (
                    (data[tid0 * nx * ny * nl + i * ny * nl + j * nl + 1] -
                     data[tid0 * nx * ny * nl + i * ny * nl + j * nl]) / ds);
            gradz[0 * nx * ny * nl + i * ny * nl + j * nl + (nl - 1)] = (
                    (data[tid0 * nx * ny * nl + i * ny * nl + j * nl + (nl - 1)] -
                     data[tid0 * nx * ny * nl + i * ny * nl + j * nl + (nl - 2)]) / ds);
            gradz[1 * nx * ny * nl + i * ny * nl + j * nl] = (
                    (data[tid1 * nx * ny * nl + i * ny * nl + j * nl + 1] -
                     data[tid1 * nx * ny * nl + i * ny * nl + j * nl]) / ds);
            gradz[1 * nx * ny * nl + i * ny * nl + j * nl + (nl - 1)] = (
                    (data[tid1 * nx * ny * nl + i * ny * nl + j * nl + (nl - 1)] -
                     data[tid1 * nx * ny * nl + i * ny * nl + j * nl + (nl - 2)]) / ds);
		}
	}
}

double lsGrid::clamp(double aa, double a1, double a2):
    return std::min(std::max(aa, a1), a2);

void lsGrid::square_vals(double* C, double* data, int nxp, int nyp, int tid){
    C[0] = data[tid * nx * ny + nxp * ny + nyp]; // C00
    C[1] = data[tid * nx * ny + nxp * ny + clamp(nyp + 1, 0, ny - 1)]; // C01
    C[2] = data[tid * nx * ny + clamp(nxp + 1, 0, nx - 1) * ny + nyp]; // C10
    C[3] = data[tid * nx * ny + clamp(nxp + 1, 0, nx - 1) * ny + clamp(nyp + 1, 0, ny - 1)];//C11
}

void lsGrid::cube_vals(double* C, double* data, int nxp, int nyp, int nlp, int tid){
    C[0] = data[tid * nx * ny * nl + nxp * ny * nl + nyp * nl + nlp]; //C000
    C[1] = data[tid * nx * ny * nl + nxp * ny * nl + clamp(nyp + 1, 0, ny - 1) * nl + nlp]; //C010
    C[2] = data[tid * nx * ny * nl + clamp(nxp + 1, 0, nx - 1) * ny * nl + nyp * nl + nlp]; //C100
    C[3] = data[tid * nx * ny * nl + clamp(nxp + 1, 0, nx - 1) * ny * nl + clamp(nyp + 1, 0, ny - 1) * nl + nlp]; //C110
    C[4] = data[tid * nx * ny * nl + nxp * ny * nl + nyp * nl + clamp(nlp + 1, 0, nl - 1)]; //C001
    C[5] = data[tid * nx * ny * nl + nxp * ny * nl + clamp(nyp + 1, 0, ny - 1) * nl + clamp(nlp + 1, 0, nl - 1)]; //C011
    C[6] = data[tid * nx * ny * nl + clamp(nxp + 1, 0, nx - 1) * ny * nl + nyp * nl + clamp(nlp + 1, 0, nl - 1)]; //C101
    C[7] = data[tid * nx * ny * nl + clamp(nxp + 1, 0, nx - 1) * ny * nl + clamp(nyp + 1, 0, ny - 1) * nl + clamp(nlp + 1, 0, nl - 1)]; //C111
}

double lsGrid::spline_interp_velo(double* U, double* Udt, double* Udx, double* Udy, double* Uds, int nxp, int nyp, int nsp0, int nsp1, double xd, double yd, double sd0, double sd1, double td, int tid1, int tid2){
    // Ux
    // --------------------------------------------------------------------
    // u0, u1, u0dt, u1dt, u0dx, u1dx, u0dy, u1dy, u0ds, u1ds
	double C[8], D[8], Ct[8], Dt[8], Cx[8], Dx[8], Cy[8], Dy[8], Cs[8], Ds[8];
    cube_vals(C, U, nxp, nyp, nsp0, tid1);
    cube_vals(D, U, nxp, nyp, nsp1, tid2);
    cube_vals(Ct, Udt, nxp, nyp, nsp0, tid1);
    cube_vals(Dt, Udt, nx, ny, ns, nxp, nyp, nsp1, tid2);
    cube_vals(Cx, Udx, nx, ny, ns, nxp, nyp, nsp0, 0);
    cube_vals(Dx, Udx, nx, ny, ns, nxp, nyp, nsp1, 1);
    cube_vals(Cy, Udy, nx, ny, ns, nxp, nyp, nsp0, 0);
    cube_vals(Dy, Udy, nx, ny, ns, nxp, nyp, nsp1, 1);
    cube_vals(Cs, Uds, nx, ny, ns, nxp, nyp, nsp0, 0);
    cube_vals(Ds, Uds, nx, ny, ns, nxp, nyp, nsp1, 1);

    // gcompute spline coefficients
    double a00 = dx * Cx[0] - (C[2] - C[0]);
    double b00 = -dx * Cx[2] + (C[2] - C[0]);
    double a10 = dx * Cx[1] - (C[3] - C[1]);
    double b10 = -dx * Cx[3] + (C[3] - C[1]);
    double a01 = dx * Cx[4] - (C[6] - C[4]);
    double b01 = -dx * Cx[6] + (C[6] - C[4]);
	double a11 = dx * Cx[5] - (C[7] - C[5]);
    double b11 = -dx * Cx[7] + (C[7] - C[5]);
    double c00 = dx * Dx[0] - (D[2] - D[0]);
    double d00 = -dx * Dx[2] + (D[2] - D[0]);
    double c10 = dx * Dx[1] - (D[3] - D[1]);
    double d10 = -dx * Dx[3] + (D[3] - D[1]);
    double c01 = dx * Dx[4] - (D[6] - D[4]);
    double d01 = -dx * Dx[6] + (D[6] - D[4]);
	double c11 = dx * Dx[5] - (D[7] - D[5]);
    double d11 = -dx * Dx[7] + (D[7] - D[5]);

    // compute new values
    double C00 = C[0] * (1. - xd) + C[2] * xd + xd * (1 - xd) * (a00 * (1 - xd) + b00 * xd);
    double C10 = C[1] * (1. - xd) + C[3] * xd + xd * (1 - xd) * (a10 * (1 - xd) + b10 * xd);
    double C01 = C[4] * (1. - xd) + C[6] * xd + xd * (1 - xd) * (a01 * (1 - xd) + b01 * xd);
    double C11 = C[5] * (1. - xd) + C[7] * xd + xd * (1 - xd) * (a11 * (1 - xd) + b11 * xd);
    double D00 = D[0] * (1. - xd) + D[2] * xd + xd * (1 - xd) * (c00 * (1 - xd) + d00 * xd);
    double D10 = D[1] * (1. - xd) + D[3] * xd + xd * (1 - xd) * (c10 * (1 - xd) + d10 * xd);
    double D01 = D[4] * (1. - xd) + D[6] * xd + xd * (1 - xd) * (c01 * (1 - xd) + d01 * xd);
    double D11 = D[5] * (1. - xd) + D[7] * xd + xd * (1 - xd) * (c11 * (1 - xd) + d11 * xd);

    // Reduce y, s and t
    double C00y = Cy[0] * (1. - xd) + Cy[2] * xd;
    double C10y = Cy[1] * (1. - xd) + Cy[3] * xd;
    double C01y = Cy[4] * (1. - xd) + Cy[6] * xd;
    double C11y = Cy[5] * (1. - xd) + Cy[7] * xd;
    double D00y = Dy[0] * (1. - xd) + Dy[2] * xd;
    double D10y = Dy[1] * (1. - xd) + Dy[3] * xd;
    double D01y = Dy[4] * (1. - xd) + Dy[6] * xd;
    double D11y = Dy[5] * (1. - xd) + Dy[7] * xd;

    double C00t = Ct[0] * (1. - xd) + Ct[2] * xd;
    double C10t = Ct[1] * (1. - xd) + Ct[3] * xd;
    double C01t = Ct[4] * (1. - xd) + Ct[6] * xd;
    double C11t = Ct[5] * (1. - xd) + Ct[7] * xd;
    double D00t = Dt[0] * (1. - xd) + Dt[2] * xd;
    double D10t = Dt[1] * (1. - xd) + Dt[3] * xd;
    double D01t = Dt[4] * (1. - xd) + Dt[6] * xd;
    double D11t = Dt[5] * (1. - xd) + Dt[7] * xd;

    double C00s = Cs[0] * (1. - xd) + Cs[2] * xd;
    double C10s = Cs[1] * (1. - xd) + Cs[3] * xd;
    double C01s = Cs[4] * (1. - xd) + Cs[6] * xd;
    double C11s = Cs[5] * (1. - xd) + Cs[7] * xd;
    double D00s = Ds[0] * (1. - xd) + Ds[2] * xd;
    double D10s = Ds[1] * (1. - xd) + Ds[3] * xd;
    double D01s = Ds[4] * (1. - xd) + Ds[6] * xd;
    double D11s = Ds[5] * (1. - xd) + Ds[7] * xd;

    // compute spline coefficients
    double a0 = dy * C00y - (C10 - C00);
    double b0 = -dy * C10y + (C10 - C00);
    double a1 = dy * C01y - (C11 - C01);
    double b1 = -dy * C11y + (C11 - C01);
    double c0 = dy * D00y - (D10 - D00);
    double d0 = -dy * D10y + (D10 - D00);
    double c1 = dy * D01y - (D11 - D01);
    double d1 = -dy * D11y + (D11 - D01);

    // compute new values
    double C0 = C00 * (1. - yd) + C10 * yd + yd * (1 - yd) * (a0 * (1 - yd) + b0 * yd);
    double C1 = C01 * (1. - yd) + C11 * yd + yd * (1 - yd) * (a1 * (1 - yd) + b1 * yd);
    double D0 = D00 * (1. - yd) + D10 * yd + yd * (1 - yd) * (c0 * (1 - yd) + d0 * yd);
    double D1 = D01 * (1. - yd) + D11 * yd + yd * (1 - yd) * (c1 * (1 - yd) + d1 * yd);

    // reduce y and t direction
    double C0s = C00s * (1. - yd) + C10s * yd;
    double C1s = C01s * (1. - yd) + C11s * yd;
    double D0s = D00s * (1. - yd) + D10s * yd;
    double D1s = D01s * (1. - yd) + D11s * yd;
    double C0t = C00t * (1. - yd) + C10t * yd;
    double C1t = C01t * (1. - yd) + C11t * yd;
    double D0t = D00t * (1. - yd) + D10t * yd;
    double D1t = D01t * (1. - yd) + D11t * yd;

    // compute spline coefficients
    double a0 = ds * C0s - (C1 - C0);
    double b0 = -ds * C1s + (C1 - C0);
    double a1 = ds * D0s - (D1 - D0);
    double b1 = -ds * D1s + (D1 - D0);

    // compute wave elevation for two time steps
    double C = C0 * (1. - sd0) + C1 * sd0 + sd0 * (1 - sd0) * (a0 * (1 - sd0) + b0 * sd0);
    double D = D0 * (1. - sd1) + D1 * sd1 + sd1 * (1 - sd1) * (a1 * (1 - sd1) + b1 * sd1);

    // compute spline coefficients
    double Ct = C0t * (1. - sd0) + C1t * sd0;
    double Dt = D0t * (1. - sd1) + D1t * sd1;
    double a = dt * Ct - (D - C);
    double b = -dt * Dt + (D - C);

    // compute final interpolated wave elevation for single point
    return C * (1. - td) + D * td + td * (1 - td) * (a * (1 - td) + b * td);
}

/* Function for trilinear interpolation on a cartesian evenly spaced mesh*/
double lsGrid::trilinear_interpolation(double* VAR0, double* VAR1, double tpt, double xpt, double ypt, double zpt) {

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
double lsGrid::trilinear_interpolation2(double* VAR0, double* VAR1, double tpt, double xpt, double ypt, double zpt) {

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
double lsGrid::bilinear_interpolation(double* VAR0, double* VAR1, double tpt, double xpt, double ypt) {

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
double lsGrid::bilinear_interpolation2(double* VAR0, double* VAR1, double tpt, double xpt, double ypt) {



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


double lsGrid::u(double tpt, double xpt, double ypt, double zpt) {
	return trilinear_interpolation2(UX0, UX1, tpt, xpt, ypt, zpt);
}

double lsGrid::v(double tpt, double xpt, double ypt, double zpt) {
	return trilinear_interpolation2(UY0, UY1, tpt, xpt, ypt, zpt);
}

double lsGrid::w(double tpt, double xpt, double ypt, double zpt) {
	return trilinear_interpolation2(UZ0, UZ1, tpt, xpt, ypt, zpt);
}

double lsGrid::eta(double tpt, double xpt, double ypt) {
	//update_bounds(xpt, ypt);
	return bilinear_interpolation2(ETA0, ETA1, tpt, xpt, ypt);
}


bool lsGrid::CheckTime(double tpt) {
	/* Checks to see if the time tpt is within the interval t0 to t1. If so, returns true*/
	if (tpt > t0 + dt) {
		std::cout << "t0: " << t0 << ", t1: " << (t0 + dt) << ", tpt: " << tpt << std::endl;
		return false;

	}
	return true;
}

// function to find if given point 
// lies inside a given rectangle or not. 
bool lsGrid::CheckBounds()
{
	if (bxmin >= domain[0] && bxmax <= domain[1] && bymin >= domain[2] && bymax <= domain[3])
		return true;
	else {
		std::cout << "Requested point outside specified grid domain. adjust the bounds of the grid and try again." << std::endl;
		exit(-1);
		return false;
	}
}

void lsGrid::update_bounds(double xpt, double ypt) {
	bxmin = std::min(xpt, bxmin);
	bxmax = std::min(xpt, bxmax);
	bymin = std::min(ypt, bymin);
	bymax = std::min(ypt, bymax);

}


//----------------------------------------------------------------------------------------------------------------------------------------


// Second order theory functions


//----------------------------------------------------------------------------------------------------------------------------------------

// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 3
void lsGrid::initialme_kinematics(Irregular& irregular) {

	dx = (domain[1] - domain[0]) / std::max(1., double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));
	//dz = (domain_end[2] - domain_start[2]) / double(NZ - 1);
	ds = 1. / std::max(1., double(nl - 1));

	double dd = omp_get_wtime();

	//omp_set_num_threads(1)
	//omp_set_num_threads(omp_get_max_threads());

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
					double zpt0max = std::max(0., zpt0 - swl);
					double zpt0min = std::min(0., zpt0); // swl already included in irregular class functions
					
					//UX0[i * ny * nl + j * nl + m] = irregular.u1(t0, xpt, ypt, zpt0min) + irregular.u2(t0, xpt, ypt, zpt0min) + PHI0_dxdz * zpt0max;
					//UY0[i * ny * nl + j * nl + m] = irregular.v1(t0, xpt, ypt, zpt0min) + irregular.v2(t0, xpt, ypt, zpt0min) + PHI0_dydz * zpt0max;
					//UZ0[i * ny * nl + j * nl + m] = irregular.w1(t0, xpt, ypt, zpt0min) + irregular.w2(t0, xpt, ypt, zpt0min) + PHI0_dzdz * zpt0max;

					std::vector<double> U2 = irregular.uvw2(t0, xpt, ypt, zpt0min);
					UX0[i * ny * nl + j * nl + m] = irregular.u1(t0, xpt, ypt, zpt0min) + U2[0] + PHI0_dxdz * zpt0max;
					UY0[i * ny * nl + j * nl + m] = irregular.v1(t0, xpt, ypt, zpt0min) + U2[1] + PHI0_dydz * zpt0max;
					UZ0[i * ny * nl + j * nl + m] = irregular.w1(t0, xpt, ypt, zpt0min) + U2[2] + PHI0_dzdz * zpt0max;
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
						double zpt1max = std::max(0., zpt1 - swl);
						double zpt1min = std::min(0., zpt1); // swl already included in irregular class functions

						//UX1[i * ny * nl + j * nl + m] = irregular.u1(t0 + dt, xpt, ypt, zpt1min) + irregular.u2(t0 + dt, xpt, ypt, zpt1min) + PHI1_dxdz * zpt1max;
						//UY1[i * ny * nl + j * nl + m] = irregular.v1(t0 + dt, xpt, ypt, zpt1min) + irregular.v2(t0 + dt, xpt, ypt, zpt1min) + PHI1_dydz * zpt1max;
						//UZ1[i * ny * nl + j * nl + m] = irregular.w1(t0 + dt, xpt, ypt, zpt1min) + irregular.w2(t0 + dt, xpt, ypt, zpt1min) + PHI1_dzdz * zpt1max;

						std::vector<double> U2 = irregular.uvw2(t0 + dt, xpt, ypt, zpt1min);
						UX1[i * ny * nl + j * nl + m] = irregular.u1(t0 + dt, xpt, ypt, zpt1min) + U2[0] + PHI1_dxdz * zpt1max;
						UY1[i * ny * nl + j * nl + m] = irregular.v1(t0 + dt, xpt, ypt, zpt1min) + U2[0] + PHI1_dydz * zpt1max;
						UZ1[i * ny * nl + j * nl + m] = irregular.w1(t0 + dt, xpt, ypt, zpt1min) + U2[0] + PHI1_dzdz * zpt1max;
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
void lsGrid::initialize_kinematics_with_ignore(Irregular& irregular) {

	dx = (domain[1] - domain[0]) / std::max(1., double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));
	//dz = (domain_end[2] - domain_start[2]) / double(NZ - 1);
	ds = 1. / std::max(1., double(nl - 1));

	double dd = omp_get_wtime();

	//omp_set_num_threads(1);
	//omp_set_num_threads(omp_get_max_threads());

#pragma omp parallel // start parallel initialization
	{
#pragma omp master
		std::cout << "Number of available threads: " << omp_get_num_threads() << std::endl;

		// Main grid
#pragma omp for collapse(2)
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
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
						double zpt0max = std::max(0., zpt0 - swl);
						double zpt0min = std::min(0., zpt0); // swl already included in irregular class functions

						//UX0[i * ny * nl + j * nl + m] = irregular.u1(t0, xpt, ypt, zpt0min) + irregular.u2(t0, xpt, ypt, zpt0min) + PHI0_dxdz * zpt0max;
						//UY0[i * ny * nl + j * nl + m] = irregular.v1(t0, xpt, ypt, zpt0min) + irregular.v2(t0, xpt, ypt, zpt0min) + PHI0_dydz * zpt0max;
						//UZ0[i * ny * nl + j * nl + m] = irregular.w1(t0, xpt, ypt, zpt0min) + irregular.w2(t0, xpt, ypt, zpt0min) + PHI0_dzdz * zpt0max;

						std::vector<double> U2 = irregular.uvw2(t0, xpt, ypt, zpt0min);
						UX0[i * ny * nl + j * nl + m] = irregular.u1(t0, xpt, ypt, zpt0min) + U2[0] + PHI0_dxdz * zpt0max;
						UY0[i * ny * nl + j * nl + m] = irregular.v1(t0, xpt, ypt, zpt0min) + U2[1] + PHI0_dydz * zpt0max;
						UZ0[i * ny * nl + j * nl + m] = irregular.w1(t0, xpt, ypt, zpt0min) + U2[2] + PHI0_dzdz * zpt0max;
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
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
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
							double zpt1max = std::max(0., zpt1 - swl);
							double zpt1min = std::min(0., zpt1); // swl already included in irregular class functions

							//UX1[i * ny * nl + j * nl + m] = irregular.u1(t0 + dt, xpt, ypt, zpt1min) + irregular.u2(t0 + dt, xpt, ypt, zpt1min) + PHI1_dxdz * zpt1max;
							//UY1[i * ny * nl + j * nl + m] = irregular.v1(t0 + dt, xpt, ypt, zpt1min) + irregular.v2(t0 + dt, xpt, ypt, zpt1min) + PHI1_dydz * zpt1max;
							//UZ1[i * ny * nl + j * nl + m] = irregular.w1(t0 + dt, xpt, ypt, zpt1min) + irregular.w2(t0 + dt, xpt, ypt, zpt1min) + PHI1_dzdz * zpt1max;

							std::vector<double> U2 = irregular.uvw2(t0 + dt, xpt, ypt, zpt1min);
							UX1[i * ny * nl + j * nl + m] = irregular.u1(t0 + dt, xpt, ypt, zpt1min) + U2[0] + PHI1_dxdz * zpt1max;
							UY1[i * ny * nl + j * nl + m] = irregular.v1(t0 + dt, xpt, ypt, zpt1min) + U2[0] + PHI1_dydz * zpt1max;
							UZ1[i * ny * nl + j * nl + m] = irregular.w1(t0 + dt, xpt, ypt, zpt1min) + U2[0] + PHI1_dzdz * zpt1max;
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
void lsGrid::initialize_surface_elevation(Irregular& irregular, double t_target) {
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
	//omp_set_num_threads(omp_get_max_threads());

#pragma omp parallel
	{
		// Main grid
#pragma omp for collapse(2)
		for (int i = 0; i < nx; i++) {

			for (int j = 0; j < ny; j++) {
				double xpt = domain[0] + dx * i;
				double ypt = domain[2] + dy * j;

				ETA0[i * ny + j] = irregular.eta1(t0, xpt, ypt) + irregular.eta2(t0, xpt, ypt) + swl;
				//ETA1[i * ny + j] = irregular.eta1(t0 + dt, xpt, ypt) + irregular.eta2(t0 + dt, xpt, ypt) + swl;

			}
		}
		if (init_only){
		  ETA1 = ETA0;
		}
		else{
#pragma omp for collapse(2)
                for (int i = 0; i < nx; i++) {

                        for (int j = 0; j < ny; j++) {
                                double xpt = domain[0] + dx * i;
                                double ypt = domain[2] + dy * j;

                                //ETA0[i * ny + j] = irregular.eta1(t0, xpt, ypt) + irregular.eta2(t0, xpt, ypt) + swl;
                                ETA1[i * ny + j] = irregular.eta1(t0 + dt, xpt, ypt) + irregular.eta2(t0 + dt, xpt, ypt) + swl;

                        }
                }
		}

	}
	dd = omp_get_wtime() - dd;

	std::cout << "Surface Elevation generated successfully. ";
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
}
// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 3
void lsGrid::initialize_surface_elevation_with_ignore(Irregular& irregular, double t_target) {

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
	//omp_set_num_threads(omp_get_max_threads());

#pragma omp parallel
	{
		// Main grid
#pragma omp for collapse(2)
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				double xpt = domain[0] + dx * i;
				double ypt = domain[2] + dy * j;
				if (!IGNORE[i * ny + j]) {
					ETA0[i * ny + j] = irregular.eta1(t0, xpt, ypt) + irregular.eta2(t0, xpt, ypt) + swl;
					//ETA1[i * ny + j] = irregular.eta1(t0 + dt, xpt, ypt) + irregular.eta2(t0 + dt, xpt, ypt) + swl;
				}
				else {
					ETA0[i * ny + j] = swl;
					//ETA1[i * ny + j] = swl;
				}
			}
		}
		if (init_only) {
		  ETA1 = ETA0;
		}
		else{
#pragma omp for collapse(2)
                for (int j = 0; j < ny; j++) {
                        for (int i = 0; i < nx; i++) {
                                double xpt = domain[0] + dx * i;
                                double ypt = domain[2] + dy * j;
                                if (!IGNORE[i * ny + j]) {
				  //ETA0[i * ny + j] = irregular.eta1(t0, xpt, ypt) + irregular.eta2(t0, xpt, ypt) + swl;
                                        ETA1[i * ny + j] = irregular.eta1(t0 + dt, xpt, ypt) + irregular.eta2(t0 + dt, xpt, ypt) + swl;
                                }
                                else {
				  //ETA0[i * ny + j] = swl;
                                        ETA1[i * ny + j] = swl;
                                }
                        }
                }
		}

	}
	dd = omp_get_wtime() - dd;

	std::cout << "Surface Elevation generated successfully. ";
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
}


// When called, updates the arrays storing surface elevation and kinematics data for timestep t0 = t1, t1 = t1+dt
void lsGrid::update(Irregular& irregular, double t_target)
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
		//omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel
		{
			// Main grid
#pragma omp for collapse(2)
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
					double xpt = domain[0] + dx * i;
					//std::cout << "processornum: " << omp_get_thread_num() << std::endl;
					double ypt = domain[2] + dy * j;
					if (!IGNORE[i * ny + j]) {
						ETA0[i * ny + j] = ETA1[i * ny + j];
						ETA1[i * ny + j] = irregular.eta1(t0 + dt, xpt, ypt) + irregular.eta2(t0 + dt, xpt, ypt) + swl;

						//double Ux1 = irregular.u1(t0 + dt, xpt, ypt, 0.0) + irregular.u2(t0 + dt, xpt, ypt, 0.0);
						//double Uy1 = irregular.v1(t0 + dt, xpt, ypt, 0.0) + irregular.v2(t0 + dt, xpt, ypt, 0.0);
						//double Uz1 = irregular.w1(t0 + dt, xpt, ypt, 0.0) + irregular.w2(t0 + dt, xpt, ypt, 0.0);

						double PHI1_dxdz = irregular.phi1_dxdz(t0 + dt, xpt, ypt);
						double PHI1_dydz = irregular.phi1_dydz(t0 + dt, xpt, ypt);
						double PHI1_dzdz = irregular.phi1_dzdz(t0 + dt, xpt, ypt);

						for (int m = 0; m < nl; m++) {
							double spt = s2tan(-1. + ds * m);
							double zpt1 = s2z(spt, ETA1[i * ny + j] - swl, water_depth);

							UX0[i * ny * nl + j * nl + m] = UX1[i * ny * nl + j * nl + m];
							UY0[i * ny * nl + j * nl + m] = UY1[i * ny * nl + j * nl + m];
							UZ0[i * ny * nl + j * nl + m] = UZ1[i * ny * nl + j * nl + m];

							double zpt1max = std::max(0., zpt1 - swl);
							double zpt1min = std::min(0., zpt1); // swl already included in irregular class functions

							//UX1[i * ny * nl + j * nl + m] = irregular.u1(t0 + dt, xpt, ypt, zpt1min) + irregular.u2(t0 + dt, xpt, ypt, zpt1min) + PHI1_dxdz * zpt1max;
							//UY1[i * ny * nl + j * nl + m] = irregular.v1(t0 + dt, xpt, ypt, zpt1min) + irregular.v2(t0 + dt, xpt, ypt, zpt1min) + PHI1_dydz * zpt1max;
							//UZ1[i * ny * nl + j * nl + m] = irregular.w1(t0 + dt, xpt, ypt, zpt1min) + irregular.w2(t0 + dt, xpt, ypt, zpt1min) + PHI1_dzdz * zpt1max;

							std::vector<double> U2 = irregular.uvw2(t0 + dt, xpt, ypt, zpt1min);
							UX1[i * ny * nl + j * nl + m] = irregular.u1(t0 + dt, xpt, ypt, zpt1min) + U2[0] + PHI1_dxdz * zpt1max;
							UY1[i * ny * nl + j * nl + m] = irregular.v1(t0 + dt, xpt, ypt, zpt1min) + U2[0] + PHI1_dydz * zpt1max;
							UZ1[i * ny * nl + j * nl + m] = irregular.w1(t0 + dt, xpt, ypt, zpt1min) + U2[0] + PHI1_dzdz * zpt1max;
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

void lsGrid::initialize_kinematics(SpectralWaveData *swd) {
	// Tell the swd object current application time...
	dx = (domain[1] - domain[0]) / std::max(1., double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));
	//dz = (domain_end[2] - domain_start[2]) / double(NZ - 1);
	ds = 1. / std::max(1., double(nl - 1));

	double dd = omp_get_wtime();

	//omp_set_num_threads(1);
	//omp_set_num_threads(omp_get_max_threads());

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
	if (init_only) {
                UX1 = UX0;
                UY1 = UY0;
                UZ1 = UZ0;
        }
	else{

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
	}
	if (dump_vtk) {
		write_vtk(false);
		update_count++;
	}
	std::cout << "Generation of domain kinematics using SWD completed. ";
	dd = omp_get_wtime() - dd;
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
	std::cout << "Interpolation can commence..." << std::endl;
}


void lsGrid::initialize_kinematics_with_ignore(SpectralWaveData* swd) {
	// Tell the swd object current application time...
	dx = (domain[1] - domain[0]) / std::max(1., double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));
	//dz = (domain_end[2] - domain_start[2]) / double(NZ - 1);
	ds = 1. / std::max(1., double(nl - 1));

	double dd = omp_get_wtime();

	//omp_set_num_threads(1);
	//omp_set_num_threads(omp_get_max_threads());

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
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
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
		 // timestep T1.
		if (init_only) {
		  UX1 = UX0;
		  UY1 = UY0;
		  UZ1 = UZ0;
		}
		else {


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
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
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

void lsGrid::initialize_surface_elevation(SpectralWaveData* swd, double t_target) {

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

	//omp_set_num_threads(omp_get_max_threads());
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
	if (init_only) {
	  ETA1 = ETA0;
	}
	else{
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
	}
	dd = omp_get_wtime() - dd;

	std::cout << "Surface Elevation generated successfully. ";
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
}

void lsGrid::initialize_surface_elevation_with_ignore(SpectralWaveData* swd, double t_target) {
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
	//omp_set_num_threads(omp_get_max_threads());

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
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
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
		if (init_only){
		  ETA1 = ETA0;
		}
		else{
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
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
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
	}
	dd = omp_get_wtime() - dd;

	std::cout << "Surface Elevation generated successfully. ";
	std::cout << "Initialization time: " << dd << " seconds." << std::endl;
}

void lsGrid::update(SpectralWaveData* swd, double t_target)
{
	// Start by checking bounds
/*
if (!disable_checkbounds){
	CheckBounds();
}*/

// new time step
	if ((t_target / dt - (t0 + 2 * dt) / dt) > 0.) {
		double new_time = dt * std::floor(t_target / dt);
		std::cout << "Time step to large. reinitializing LSgrid." << std::endl;
		if (ignore_domain){
			initialize_surface_elevation_with_ignore(swd, new_time);
			initialize_kinematics_with_ignore(swd);
		}
		else{
			initialize_surface_elevation(swd, new_time);
			initialize_kinematics(swd);
		}
	}
	else {
		t0 += dt;
		// Updating surface elevations
		double dd = omp_get_wtime();
		//omp_set_num_threads(1);
		//omp_set_num_threads(omp_get_max_threads());

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
// switch order of i and j on purpose since this works much betten when ignore is on and wave propagate from the x boundary only
#pragma omp for collapse(2)
			for (int j = 0; j < ny; j++) { 
				for (int i = 0; i < nx; i++) {
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
void lsGrid::write_vtk(bool endtime) {
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
void lsGrid::export_vtu(FILE* fp, bool last)
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
void lsGrid::set_ignore()
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
