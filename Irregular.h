#ifndef Irregular_H
#define Irregular_H
class Irregular {
private:
	// Variables
	int nfreq, ndir, extmet, bandwidth;
	int pertmet = 0;
	double ampl, depth, s, mtheta, tofmax, fpoint[2];
	
	// Declaration of pointers where data will be stored
	double* omega;
	double* Ampspec;
	double* k;
	double* thetaA;
	double* D;
	double* phase;
	double* dsum2;

public:
	// First order
	double eta1(double t, double xx, double yy); // wave elevation
	double u1(double t, double xx, double yy, double zz); // velocity component x
	double v1(double t, double xx, double yy, double zz); // velocity component y
	double w1(double t, double xx, double yy, double zz); // velocity component z
	double dp1(double t, double xx, double yy, double zz); // dynamic pressure component	

	// second order
	double eta2(double t, double xx, double yy);
	double u2(double t, double xx, double yy, double zz);
	double v2(double t, double xx, double yy, double zz);
	double w2(double t, double xx, double yy, double zz);

	// gradients
	double phi1_dxdz(double t, double xx, double yy); /* vertical velocity gradient at z=0 for velocity component U */
	double phi1_dydz(double t, double xx, double yy); /* vertical velocity gradient at z=0 for velocity component V */
	double phi1_dzdz(double t, double xx, double yy); /* vertical velocity gradient at z=0 for velocity component W */


	// main kinematics functions
	double eta(double t, double x, double y);
	double u(double t, double x, double y, double z); 
	double v(double t, double x, double y, double z);
	double w(double t, double x, double y, double z);


};

#endif