
// Variables
int nfreq, ndir, wavetype, extmet, pertmet, meth, bandwidth, n_timesteps, rampswitch, normalizeA, spreadfunc;
int NX, NY, NZ, NXL, NYL, NZL;
double ampl, depth, s, mtheta, tofmax, fpoint[2], trampdata[3], xrampdata[3], yrampdata[3];
double fp, alpha_z, alpha_u, ramp_time, x_pos, y_pos, current_speed, wave_length, wave_height;


#define PI 3.1415926535897
#define G 9.81
#define RHO 1025.0
//fftw_plan p;

// Declaration of pointers where data will be stored
double* w;
double* Ampspec;
double* k;
double* thetaA;
double* D;
double* phas;
double* dsum2;
double* PD_time;
double* PD_ampl;
double* PD_velo;
double* PD_eta;
double* domainsize;
double* index;

// Functions for reading and processing waveinput.dat
int read_inputdata_v2();
int read_inputdata();