#ifndef Wavemaker_H
#define Wavemaker_H

class Wavemaker {

private:
	int findNearestNeighbourIndex(double value, double* x, int len);
	double interp1(double* x, int x_tam, double* y, double xx);
public:
	bool initialized = false;
	double* PD_time;
	double* PD_ampl;
	double* PD_velo;
	double* PD_eta;
	double u_piston(double t);
	double wave_elev_piston(double t);

	double alpha_z, alpha_u;
	int n_timesteps;

};


#endif
