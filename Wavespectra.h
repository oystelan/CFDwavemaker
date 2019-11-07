#ifndef Wavespectra_H
#define Wavespectra_H


/* This class is still under contruction*/

class Wavespectra {
private:
	double gamma(double f, double Gamma, bool gg);
	void jonswap3(double* S, double* omega, int nfreq, double hs, double tp, double gamma);
	void torsethaugen2004(double* S, double* omega, int nfreq, double hs, double tp);
	void torsethaugen1996(double* S, double* omega, double hs, double tp);
	void PM(double* S, double* omega, int nfreq, double hs, double tp);
	void spreading_uniform(double* D, int nfreq, int ndir);
	void spreading_cos_theta_n(double* D, int nfreq, int ndir, double s, double mtheta, double* theta);
	void spreading_cos_theta05_2s(double* D, int nfreq, int ndir, double s, double mtheta, double* theta);
	void spreading_ewans(double tp, double mtheta);

public:
	double* S;
	double* D;
	double* omega;
	double* theta;
	double* k;
	int nfreq;
	int ndir;

	//todo: implement functions generates an amplitude spectrum, populating the following arrays: omega, theta, phase and A. 

};
#endif
