#ifndef Wavespectra_H
#define Wavespectra_H

#include "Irregular.h"
#include <vector>


// structure for storage of input parameters required to initialize an event
struct WaveSpecData {
	// frequency def
	double w_min;
	double w_max;
	double dw;
	// random seed
	double random_seed;
	// wave energy spectrum def
	std::string spectrum_type;
	double hs;
	double tp;
	double gamma;
	//spreading function
	std::string spreading_type;
	double spreadfactor;
	double sampling_mode; // either integration or single_sample
	double theta_min;
	double theta_max;
	double dtheta;
};


/* This class is still under contruction*/

class Wavespectra {
private:
	double gamma_spectrum(double f, double Gamma, bool gg);
	void jonswap3(double* S, double* omega, int nfreq, double hs, double tp, double gamma);
	void torsethaugen2004(double* S, double* omega, int nfreq, double hs, double tp);
	void torsethaugen1996(double* S, double* omega, double hs, double tp);
	void PM(double* S, double* omega, int nfreq, double hs, double tp);
	void spreading_uniform(std::vector<double>& D, int nfreq, int ndir);
	void spreading_cos_theta_n(std::vector<double>& D, int nfreq, int ndir, double s, double mtheta, std::vector<double>& theta);
	void spreading_cos_theta05_2s(std::vector<double>& D, int nfreq, int ndir, double s, double mtheta, std::vector<double>& theta);
	void spreading_ewans(std::vector<double>& D, int nfreq, int ndir, double tp, double mtheta, std::vector<double>& theta);
	void randomize_amplitudes(double* ampl, int nampl);


public:
	std::vector<double> omega;
	std::vector<double> A;
	std::vector<double> k;
	std::vector<double> theta;
	std::vector<double> phase;

	// input variables required for initialization
	std::string spectrum_type = "jonswap";
	std::string spreading_type = "cosn";
	int nfreq;
	int ndir;
	double Hs, Tp, s, n, rand_seed;
	
	void initialize_seastate(Irregular& irregular, WaveSpecData& specdata);



	//todo: implement functions generates an amplitude spectrum, populating the following arrays: omega, theta, phase and A. 

};
#endif
