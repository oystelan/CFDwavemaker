#include "Wavemaker.h"
#include <cfloat>

// Linear interpolation function

int Wavemaker::findNearestNeighbourIndex(double value, double* x, int len)
{
	double dist;
	int idx;
	int i;
	idx = -1;
	dist = DBL_MAX;
	for (i = 0; i < len; i++) {
		double newDist = value - x[i];
		if (newDist >= 0 && newDist < dist) {
			dist = newDist;
			idx = i;
		}
	}
	return idx;
}
double Wavemaker::interp1(double* x, int x_tam, double* y, double xx)
{
	double dx, dy, slope, intercept, yy;
	int indiceEnVector;


	indiceEnVector = findNearestNeighbourIndex(xx, x, x_tam);
	if (indiceEnVector != -1) {
		dx = x[indiceEnVector + 1] - x[indiceEnVector];
		dy = y[indiceEnVector + 1] - y[indiceEnVector];
		slope = dy / dx;
		intercept = y[indiceEnVector] - x[indiceEnVector] * slope;
		yy = slope * xx + intercept;
	}
	else
		yy = DBL_MAX;

	return yy;
}


/* Horizontal velocity taken directly from the timeseries*/
double Wavemaker::u_piston(double t) {
	double ux = interp1(PD_time, n_timesteps, PD_velo, t);

	return ux + alpha_u * ux;
}

/* Wave elevation taken directly from piston timeseries*/
double Wavemaker::wave_elev_piston(double t) {
	return alpha_z * interp1(PD_time, n_timesteps, PD_eta, t);
}