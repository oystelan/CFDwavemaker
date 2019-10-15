#ifndef Irregular_H
#define Irregular_H
class Irregular {
private:
	// Variables
	int nfreq, ndir, extmet, pertmet, bandwidth;
	double ampl, depth, s, mtheta, tofmax, fpoint[2];
	
	// Declaration of pointers where data will be stored
	double* omega;
	double* Ampspec;
	double* k;
	double* thetaA;
	double* D;
	double* phi;
	double* dsum2;


	// First order
	double eta(double t, double xx, double yy); // wave elevation
	double u(double t, double xx, double yy, double zz); // velocity component x
	double v(double t, double xx, double yy, double zz); // velocity component y
	double w(double t, double xx, double yy, double zz); // velocity component z
	double dp(double t, double xx, double yy, double zz); // dynamic pressure component
	
	// gradients
	double phi_dxdz(double t, double xx, double yy); /* vertical velocity gradient at z=0 for velocity component U */
	double phi_dydz(double t, double xx, double yy); /* vertical velocity gradient at z=0 for velocity component V */
	double phi_dzdz(double t, double xx, double yy); /* vertical velocity gradient at z=0 for velocity component W */

	// second order
	double eta2(double t, double xx, double yy);
	double u2(double t, double xx, double yy, double zz);
	double v2(double t, double xx, double yy, double zz);
	double w2(double t, double xx, double yy, double zz);

public:


};

#endif