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

void Wavemaker::gradient1D(double* time, double* data, double* graddt){
    // function for computation of temporal gradients using central difference
    
	#pragma omp for
	for (int i = 1; i < (n_timesteps-1); i++) {
		double dt = time[i+1]-time[i-1];
		graddt[i] = (data[i+1] - data[i-1]) / (2*dt);
	}
	// then we add start and end
	double dt2 = time[1]-time[0];
	graddt[0] = (data[1]-data[0])/dt2;
	dt2 = time[n_timesteps-1]-time[n_timesteps-2];
	graddt[n_timesteps-1] = (data[n_timesteps-1]-data[n_timesteps-2])/dt2;
    
}

double Wavemaker::interp1_spline(double* x, int nx, double* y, double* grady, double xx)
{
	double dx, yy, xd;
	int indiceEnVector;


	indiceEnVector = findNearestNeighbourIndex(xx, x, nx);
	if (indiceEnVector != -1) {
		xd = (xx-x[indiceEnVector])/(x[indiceEnVector+1]-x[indiceEnVector]);
		dx = x[indiceEnVector + 1] - x[indiceEnVector];
		double a = dx * grady[indiceEnVector] - (y[indiceEnVector+1] - y[indiceEnVector]);
    	double b = -dx * grady[indiceEnVector+1] + (y[indiceEnVector+1] - y[indiceEnVector]);
		
		yy = y[indiceEnVector] * (1. - xd) + y[indiceEnVector+1] * xd + xd * (1 - xd) * (a * (1 - xd) + b * xd);
	}
	else
		yy = DBL_MAX;

	return yy;
}


/* Horizontal velocity taken directly from the timeseries*/
double Wavemaker::u_piston(double t) {
	double ux = interp1(PD_time, n_timesteps, PD_velo, t);
	//double ux = interp1_spline(PD_time, n_timesteps, PD_velo, PD_velo_grad, t);
	return ux;
}

/* Wave elevation taken directly from piston timeseries*/
double Wavemaker::wave_elev_piston(double t) {
	return interp1(PD_time, n_timesteps, PD_eta, t);
}

/* Wave maker position, interpolating timeseries*/
double Wavemaker::position_piston(double t) {
	return interp1(PD_time, n_timesteps, PD_ampl, t);
}