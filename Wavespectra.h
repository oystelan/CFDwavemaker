#ifndef Wavespectra_H
#define Wavespectra_H


class Wavespectra {
private:
	double* D;
	double gamma(double f, double Gamma, bool gg);
	void jonswap3(double* S, double* omega, int nfreq, double hs, double tp, double gamma);
	void torsethaugen2004(double* S, double* omega, int nfreq, double hs, double tp);
	void torsethaugen1996(double* S, double* omega, double hs, double tp);
	void PM(double* S, double* omega, int nfreq, double hs, double tp);
	void spreading_uniform(double* D, int nfreq, int ndir);
	void spreading_cos_theta_n(double* D, int nfreq, int ndir, double s, double mtheta, double* theta);
	void spreading_cos_theta05_2s(double* D, int nfreq, int ndir, double s, double mtheta, double* theta);
	void spreading_ewans(double* D, int nfreq, int ndir, double mtheta, double* theta, double* omega, double tp);
};
#endif
