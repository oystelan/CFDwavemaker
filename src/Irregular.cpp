#include "Irregular.h"
#include <math.h>
#include <algorithm>
#include <numeric>
#include <omp.h>



double Irregular::eta(double t, double x, double y) {
	switch(order){
		case 2:
			//return eta2(t,x,y);
			return (swl + eta1(t, x, y) + eta2(t, x, y));
		case 3:
			//return eta2(t,x,y);
			return (swl + eta1(t, x, y) + eta2_boost(t, x, y));
		default:
			//std::cout << eta1(t, x, y) << std::endl;
			return (swl + eta1(t, x, y));
			
	}
};

double Irregular::u(double t, double x, double y, double z) {
	if (vertical_theory == 1) return u_smit(t, x, y, z);

	switch (order) {
		case 3:
			{double z0 = std::min(0., z - swl);
			//return u2(t,x,y,z0);}
			//std::cout << "ape" << std::endl;
			return (u1(t, x, y, z0) + u2_boost(t, x, y, z0) + phi1_dxdz(t, x, y) * std::max(0., z - swl));}
			//std::vector<double> U = uvw2(t, x, y, z);
			//return U[0];} 
		case 2:
			{double z0 = std::min(0., z - swl);
			//return u2(t,x,y,z0);}
			//std::cout << "ape" << std::endl;
			return (u1(t, x, y, z0) + u2(t, x, y, z0) + phi1_dxdz(t, x, y) * std::max(0., z - swl));}
			//std::vector<double> U = uvw2(t, x, y, z);
			//return U[0];} 
		case 1:
			{z = std::min(0., z - swl);
			//std::cout << "ape2" << std::endl;
			return (u1(t, x, y, z) + u2(t, x, y, z));} 
			
		default:
			{return u1(t, x, y, z - swl);} 
			
	}
};
double Irregular::v(double t, double x, double y, double z) {
	if (vertical_theory == 1) return v_smit(t, x, y, z);
	switch (order) {
		case 3:
			{double z0 = std::min(0., z - swl);
			//return v2(t,x,y,z0);}
			return (v1(t, x, y, z0) + v2_boost(t, x, y, z0) + phi1_dydz(t, x, y) * std::max(0., z - swl)); }
			//std::vector<double> U = uvw2(t, x, y, z);
			//return U[1];} 
		case 2:
			{double z0 = std::min(0., z - swl);
			//return v2(t,x,y,z0);}
			return (v1(t, x, y, z0) + v2(t, x, y, z0) + phi1_dydz(t, x, y) * std::max(0., z - swl)); }
			//std::vector<double> U = uvw2(t, x, y, z);
			//return U[1];}  
		case 1:
			z = std::min(0., z - swl);
			return (v1(t, x, y, z) + v2(t, x, y, z));
			
		default:
			return v1(t, x, y, z - swl);
			
	}
};
double Irregular::w(double t, double x, double y, double z) {
	if (vertical_theory == 1) return w_smit(t, x, y, z);
	switch (order) {
		case 3:
			{double z0 = std::min(0., z - swl);
			return (w1(t, x, y, z0) + w2_boost(t, x, y, z0) + phi1_dzdz(t, x, y) * std::max(0., z - swl)); }
			//return w2(t,x,y,z0);}  
		case 2:
			{double z0 = std::min(0., z - swl);
			return (w1(t, x, y, z0) + w2(t, x, y, z0) + phi1_dzdz(t, x, y) * std::max(0., z - swl)); }
			//return w2(t,x,y,z0);}   
		case 1:
			z = std::min(0., z - swl);
			return (w1(t, x, y, z) + w2(t, x, y, z));
			
		default:
			return w1(t, x, y, z - swl);
			
	}
};

/* velocity profile for horizontal component. For implmentation of sloping bottom*/
double Irregular::profileX(int ind, double x, double y, double z) {
	switch (sloping_bottom) {
	case 0:
		return (cosh(k[ind] * (z + depth)) / cosh(k[ind] * depth));
	case 1:
		return 0.0;
	}
	return 0.0;
}

double Irregular::profileZ(int ind, double x, double y, double z) {
	switch (sloping_bottom) {
	case 0:
		return (sinh(k[ind] * (z + depth)) / cosh(k[ind] * depth));
	case 1:
		return 0.0;
	}
	return 0.0;
}

/* Smit et al. (2017) s-coordinate stretched profiles.
   The vertical argument k*(z+h) is replaced by k*h*(h+z)/(h+eta),
   which maps the instantaneous water column [−h, eta] onto the
   reference column [−h, 0]. This is equivalent to Wheeler stretching. */
double Irregular::profileX_s(int ind, double z, double eta) {
	double sf = depth * (depth + z) / (depth + eta);
	return cosh(k[ind] * sf) / cosh(k[ind] * depth);
}

double Irregular::profileZ_s(int ind, double z, double eta) {
	double sf = depth * (depth + z) / (depth + eta);
	return sinh(k[ind] * sf) / cosh(k[ind] * depth);
}

double Irregular::profileX_s2(double k_nm, double z, double eta) {
	double sf = depth * (depth + z) / (depth + eta);
	return cosh(k_nm * sf) / cosh(k_nm * depth);
}

double Irregular::profileZ_s2(double k_nm, double z, double eta) {
	double sf = depth * (depth + z) / (depth + eta);
	return sinh(k_nm * sf) / cosh(k_nm * depth);
}

double Irregular::dp(double t, double x, double y, double z) {
	if (vertical_theory == 1) return dp_smit(t, x, y, z);
	switch (order) {
		case 2:
			{double z0 = std::min(0., z - swl);
			return (dp1(t, x, y, z0) + dp2(t, x, y, z0) + dp1_dz(t, x, y) * std::max(0., z - swl)); }
			//return w2(t,x,y,z0);}   
		case 1:
			z = std::min(0., z - swl);
			return (dp1(t, x, y, z) + dp2(t, x, y, z));
			
		default:
			return dp1(t, x, y, z - swl);
			
	}
};

/* Function for printing wave components (sometimes useful for QA)*/
void Irregular::print()
{
	std::cout << "Irregular Waves" << std::endl;
	std::cout << "===============" << std::endl;
	std::cout << "\n";
	std::cout << "Number of frequencies: " << nfreq << std::endl;
	std::cout << "Number of wave directions " << ndir << std::endl;
	std::cout << "\n";
	std::cout << "Wave components:\n";
	std::cout << "----------------\n";
	std::cout << "|  Omega (rad/s)   |  Amplitude (m)  |  Wave number k  |   phase (rad)   | \n";
	for (int i = 0; i < (nfreq * ndir); i++) {
		std::cout << omega[i] << " " << A[i] << " " << k[i] << " " << phase[i] << std::endl;
	}
	std::cout << std::endl;
}

/* First order wave elevation */
double Irregular::eta1(double t, double xx, double yy) {

	double welev = 0.0;
	double phi;
	//std::cout << "ndir" << ndir << std::endl;
	//std::cout << "nfreq" << nfreq << std::endl;
	//std::cout << "mtheta" << mtheta << std::endl;
	
	for (int i = 0; i < ndir * nfreq; i++) {

		phi = omega[i] * tofmax + phase[i];
		welev += A[i] * cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) 
			+ sin(theta[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi) *
			freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
		//std::cout << i << ": A" << A[i] << std::endl;
	}


	return welev;

}

double Irregular::freq_ramp(double t, double ramptime, double tstart, double duration) {
	double rampdown = std::max(std::min(1. - ((t - tstart- duration + ramptime) / (ramptime + SMALL)),1.), 0.);
	double rampup = std::max(std::min(((t - tstart) / (ramptime + SMALL)), 1.), 0.);
	return rampup*rampdown;
}

int Irregular::bw_limiter(int i1) {
	int cc = i1;
	while ((omega[cc]-omega[i1] < dw_bandwidth) && (cc < nfreq * ndir)) {	
		cc++;
	}
	return cc;
}

void Irregular::calculate_bwindices() {
	bwlim.clear();
	for (int i = 0; i < nfreq * ndir ; i++) {		
		bwlim.push_back(bw_limiter(i));
		//std::cout << "i: " << i << ", bwlim: " << bw_limiter(i) << ", omega: " << omega[i] << ", omegalim: " << omega[bw_limiter(i)-1] << std::endl;
	}
}


double Irregular::eta2(double t, double xx, double yy)
{
    const int N = nfreq * ndir;
    const double dxp = xx - fpoint[0];
    const double dyp = yy - fpoint[1];

    if ((int)cth_cache.size() != N) {
        rebuild_second_order_cache();
    }

    double eta2_total = 0.0;

    std::vector<double> cphi(N), sphi(N), ramp(N);

    for (int i = 0; i < N; ++i) {
        const double ph =
            k[i] * (cth_cache[i] * dxp + sth_cache[i] * dyp)
            - omega[i] * t
            + phase0_cache[i];

        cphi[i] = std::cos(ph);
        sphi[i] = std::sin(ph);
        ramp[i] = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
    }

    // diagonal terms
    for (int i = 0; i < N; ++i) {
        if (omega[i] > dw_cutoff) {
            continue;
        }

        const double ki = k[i];
        const double Rni = R_cache[i];
        const double Rsni = Rs_cache[i];

        const double D_nn_plus =
            (6.0 * Rsni * Rsni * (ki * ki - Rni * Rni)) /
            (2.0 * Rsni * Rsni - ki * std::tanh(2.0 * ki * depth));

        const double alpha_nn_plus =
            (D_nn_plus - (ki * ki - Rni * Rni)) / Rni + 2.0 * Rni;

        if (2.0 * ki <= kmax_cache) {
            const double cos_2phi_i = cphi[i] * cphi[i] - sphi[i] * sphi[i];
            eta2_total += 0.25 * A[i] * A[i] * alpha_nn_plus * cos_2phi_i * ramp[i];
        }
    }

    // off-diagonal terms from pair cache
    for (const auto& pair : pair_cache) {
        const int i = pair.i;
        const int m = pair.m;

        const double amp = 0.5 * A[i] * A[m];
        const double ri = ramp[i];  // preserves your original weighting convention

        if (pair.has_plus) {
            const double cos_phi_plus =
                cphi[i] * cphi[m] - sphi[i] * sphi[m];
            eta2_total += amp * pair.alpha_nm_plus * cos_phi_plus * ri;
        }

        if (pair.has_minus) {
            const double cos_phi_minus =
                cphi[i] * cphi[m] + sphi[i] * sphi[m];
            eta2_total += amp * pair.alpha_nm_minus * cos_phi_minus * ri;
        }
    }

    return eta2_total;
}

/* Function which precomputes all second order terms which are not a function of space and time. Speeds up second order calculation substantially */
void Irregular::second_order_init() {
	// we start by allocating memory
	eta2_diag.resize(nfreq*ndir);
	eta2_sum.resize(nfreq*ndir*nfreq*ndir);
	eta2_diff.resize(nfreq*ndir*nfreq*ndir);
	uvw2_diag.resize(nfreq*ndir);
	uvw2_sum.resize(nfreq*ndir*nfreq*ndir);
	uvw2_diff.resize(nfreq*ndir*nfreq*ndir);
	double kmax = *max_element(k.begin(), k.end());
	
	for (int i = 0; i < nfreq*ndir; i++) {
		//double eta2_t = 0;
		if (omega[i] <= dw_cutoff) {
			double thi = theta[i] + (mtheta * PI / 180.);
			//double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];
			double Rn = k[i] * tanh(k[i] * depth);
			double Rsn = sqrt(Rn);

			// Diagonal terms
			double D_nn_plus = (6. * Rsn * Rsn * (k[i] * k[i] - Rn * Rn)) / (2. * Rsn * Rsn - k[i] * tanh(2. * k[i] * depth));
			double alpha_nn_plus = (D_nn_plus - (k[i] * k[i] - Rn * Rn)) / Rn + 2. * Rn;

			if (2 * k[i] <= kmax) {
				eta2_diag[i] = 0.25 * A[i] * A[i] * alpha_nn_plus;// * cos(2 * phi_i);
				uvw2_diag[i] = 0.25*(A[i] * g/ omega[i])*(A[i] * g / omega[i])*(D_nn_plus/(2*omega[i]))*2*k[i];
			}
			



			// Off-diagonal (sum over half only due to symmetry)
			for (int m = i + 1; m < nfreq * ndir; m++) {
				
				if ((omega[m] <= dw_cutoff) && abs(omega[i] - omega[m]) <= dw_bandwidth) {
					double thm = theta[m] + (mtheta * PI / 180.);
					//double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];

					double gamma_nm = cos(thi - thm);
					double k_nm_plus = sqrt(k[i] * k[i] + k[m] * k[m] + (2. * k[i] * k[m] * gamma_nm));
					double k_nm_minus = sqrt(k[i] * k[i] + k[m] * k[m] - (2. * k[i] * k[m] * gamma_nm));

					double Rm = k[m] * tanh(k[m] * depth);
					double Rsm = sqrt(Rm);
					double Rsnm2plus = pow(Rsn + Rsm, 2.);
					double Rsnm2minus = pow(Rsn - Rsm, 2.);


					if ((k[i] + k[m]) <= kmax) {
						double D_nm_plus = (2. * Rsnm2plus * (k[i] * k[m] * gamma_nm - Rn * Rm)) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth)) +
							((Rsn + Rsm) * (Rsn * (k[m] * k[m] - Rm * Rm) + Rsm * (k[i] * k[i] - Rn * Rn))) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth));
						double alpha_nm_plus = (D_nm_plus - (k[i] * k[m] * gamma_nm - Rn * Rm)) / sqrt(Rn * Rm) + (Rn + Rm);
						eta2_sum[i*(nfreq * ndir)+ m] = 0.5 * A[i] * A[m] * alpha_nm_plus;// * cos(phi_i + phi_m);
						uvw2_sum[i*(nfreq * ndir)+ m] = 0.5*(A[i] * g / omega[i]) * (A[m] * g / omega[m]) * (D_nm_plus/(omega[i]+omega[m]));
						
					}
					
					if (k[i] != k[m]) {
						double D_nm_minus = (2. * Rsnm2minus * (k[i] * k[m] * gamma_nm + Rn * Rm)) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth)) 
							+ ((Rsn - Rsm) * (-Rsn * (k[m] * k[m] - Rm * Rm) + Rsm * (k[i] * k[i] - Rn * Rn))) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth));
						double alpha_nm_minus = (D_nm_minus - (k[i] * k[m] * gamma_nm + Rn * Rm)) / sqrt(Rn * Rm) + (Rn + Rm);
						eta2_diff[i*(nfreq * ndir)+ m] = 0.5 * A[i] * A[m] * alpha_nm_minus;// * cos(phi_i - phi_m); // difference term
						uvw2_diff[i*(nfreq * ndir)+ m] = 0.5*(A[i] * g / omega[i]) * (A[m] * g / omega[m]) * (D_nm_minus/(omega[i]-omega[m]));
						
					}
					
				}

			}
		}
	}
}

/* Second order wave elevation */
double Irregular::eta2_boost(double t, double xx, double yy) {
	double eta2_ti = 0.;
	#pragma omp parallel for reduction(+:eta2_ti)
	for (int i = 0; i < nfreq*ndir; i++) {
		double eta2_t = 0;
		
		double thi = theta[i] + (mtheta * PI / 180.);
		double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];
		
		eta2_t += eta2_diag[i]*cos(2 * phi_i);

		// Off-diagonal (sum over half only due to symmetry)
		for (int m = i + 1; m < nfreq * ndir; m++) {
			double thm = theta[m] + (mtheta * PI / 180.);
			double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];

			eta2_t += eta2_sum[i*(nfreq * ndir)+ m] * cos(phi_i + phi_m); // sum term			
			eta2_t += eta2_diff[i*(nfreq * ndir)+ m] * cos(phi_i - phi_m); // difference term
		}
		eta2_ti += eta2_t*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return eta2_ti;
}


double Irregular::u1(double t, double xx, double yy, double zz) {

	double usum = 0.0;
	double phi;
	//(cosh(k[i] * (zz + depth)) / sinh(k[i] * depth))
	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		usum +=  A[i] * g * k[i] / omega[i] * profileX(i, xx, yy, zz)
			* cos(theta[i] + (mtheta * PI / 180.)) * cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.))
				* (xx - fpoint[0]) + sin(theta[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi)*
			freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}

	return usum;

}


/* horizontal velocity U for a sinus Irregular::omegaave */
double Irregular::v1(double t, double xx, double yy, double zz) {

	double vsum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		vsum += A[i] * g * k[i] / omega[i] * profileX(i, xx, yy, zz)
			* sin(theta[i] + (mtheta * PI / 180.)) * cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0])
				+ sin(theta[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi)*
			freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}

	return vsum;

}


/* vertical velocity for a sinus wave */
double Irregular::w1(double t, double xx, double yy, double zz) {

	double wsum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		wsum += A[i] * g * k[i] / omega[i] * profileZ(i, xx, yy, zz)
			* sin(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[i] + (mtheta 
				* PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi)*
			freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}

	return wsum;

}

/* Second order horizontal velocity U for a sinus wave */
double Irregular::u2(double t, double xx, double yy, double zz) {

	double usum2_t = 0;
	
	for (int i = 0; i < nfreq * ndir; i++) {
		double usum2 = 0;
		double thi = theta[i]+ (mtheta * PI / 180.);
		double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];
		double Rn = k[i] * tanh(k[i] * depth);
		double Rsn = sqrt(Rn);

		// Diagonal terms
		double D_nn_plus = ( 6.*Rsn*Rsn * (k[i] * k[i] - Rn * Rn))/(2.* Rsn* Rsn - k[i] * tanh(2*k[i] * depth));		
		double AA = 0.25*(A[i] * g/ omega[i])*(A[i] * g / omega[i])*(D_nn_plus/(2*omega[i]));
		usum2 += AA*(2*k[i]*cos(thi)) * cos(2*phi_i) * cosh(2*k[i]*(zz+depth)) / cosh(2*k[i]*depth);

		// Off-diagonals
		for (int m = i + 1; m < nfreq * ndir; m++) {
			double thm = theta[m]+ (mtheta * PI / 180.);
			
			double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];

			double gamma_nm = cos(thi - thm);
			double k_nm_plus = sqrt(k[i] * k[i] + k[m] * k[m] + (2. * k[i] * k[m] * gamma_nm));
			double k_nm_minus = sqrt(k[i] * k[i] + k[m] * k[m] - (2. * k[i] * k[m] * gamma_nm));

			double Rm = k[m] * tanh(k[m] * depth);

			double Rsm = sqrt(Rm);
			
			double Rsnm2plus = pow(Rsn + Rsm, 2.);
			double Rsnm2minus = pow(Rsn - Rsm, 2.);

			double D_nm_plus = (2. * Rsnm2plus * (k[i] * k[m] * gamma_nm - Rn * Rm)) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth)) + 
			((Rsn + Rsm) * (Rsn*(k[m]*k[m] - Rm * Rm) + Rsm*(k[i] * k[i] - Rn * Rn))) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth));

			double D_nm_minus = (2. * Rsnm2minus * (k[i] * k[m] * gamma_nm + Rn * Rm)) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth)) + ((Rsn - Rsm) * (-Rsn*(k[m]*k[m] - Rm * Rm) + Rsm*(k[i]*k[i] - Rn * Rn))) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth));

			double AA = 0.5*(A[i] * g / omega[i]) * (A[m] * g / omega[m]);
			double BB = (D_nm_plus/(omega[i]+omega[m]));
			double coshcosh_plus = cosh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth);

			// Sum terms
			usum2 += AA * BB * (k[i]*cos(thi) + k[m]*cos(thm)) * cos(phi_i + phi_m) * coshcosh_plus;
	
			if (k[i] != k[m]){
				double CC = (D_nm_minus/(omega[i]-omega[m])); 
				double coshcosh_minus = cosh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth);
				
				// Diff terms
				usum2 += AA * CC * (k[i]*cos(thi)-k[m]*cos(thm)) * cos(phi_i-phi_m) * coshcosh_minus;
			}	
		}
		usum2_t += usum2*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return usum2_t;
}

/* Second order horizontal velocity U for a sinus wave */
double Irregular::u2_boost(double t, double xx, double yy, double zz) {

	double usum2_t = 0;
	#pragma omp parallel for reduction(+:usum2_t)
	for (int i = 0; i < nfreq * ndir; i++) {
		double usum2 = 0;
		double thi = theta[i]+ (mtheta * PI / 180.);
		double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];

		// Diagonal terms
		usum2 += uvw2_diag[i] * cos(thi) * cos(2*phi_i) * cosh(2*k[i]*(zz+depth)) / cosh(2*k[i]*depth);

		// Off-diagonals
		for (int m = i + 1; m < nfreq * ndir; m++) {
			double thm = theta[m]+ (mtheta * PI / 180.);
			double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];
			double gamma_nm = cos(thi - thm);
			double k_nm_plus = sqrt(k[i] * k[i] + k[m] * k[m] + (2. * k[i] * k[m] * gamma_nm));
			double k_nm_minus = sqrt(k[i] * k[i] + k[m] * k[m] - (2. * k[i] * k[m] * gamma_nm));
			double coshcosh_plus = cosh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth);
			double coshcosh_minus = cosh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth);

			// Sum terms
			usum2 +=  uvw2_sum[i*(nfreq * ndir)+ m] * (k[i]*cos(thi) + k[m]*cos(thm)) * cos(phi_i + phi_m) * coshcosh_plus;
	
			// Diff terms
			usum2 += uvw2_diff[i*(nfreq * ndir)+ m] * (k[i]*cos(thi) - k[m]*cos(thm)) * cos(phi_i - phi_m) * coshcosh_minus;
		}
		usum2_t += usum2*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return usum2_t;
}


/* Second order horizontal velocity V for a sinus wave */
double Irregular::v2(double t, double xx, double yy, double zz) {

	//double eta1_t = 0;
	double vsum2_t = 0;

	for (int i = 0; i < nfreq * ndir; i++) {
		double vsum2 = 0;		
		double thi = theta[i]+ (mtheta * PI / 180.);
		double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];
		double Rn = k[i] * tanh(k[i] * depth);
		double Rsn = sqrt(Rn);

		// Diagonal terms
		double D_nn_plus = ( 6.*Rsn*Rsn * (k[i] * k[i] - Rn * Rn))/(2.* Rsn* Rsn - k[i] * tanh(2*k[i] * depth));		
		double AA = 0.25*(A[i] * g/ omega[i])*(A[i] * g / omega[i])*(D_nn_plus/(2*omega[i]));
		vsum2 += AA*(2*k[i]*sin(thi)) * cos(2*phi_i) * cosh(2*k[i]*(zz+depth)) / cosh(2*k[i]*depth); 

		for (int m = i + 1; m < nfreq * ndir; m++) {
			double thm = theta[m]+ (mtheta * PI / 180.);
			
			double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];

			double gamma_nm = cos(thi - thm);
			double k_nm_plus = sqrt(k[i] * k[i] + k[m] * k[m] + (2. * k[i] * k[m] * gamma_nm));
			double k_nm_minus = sqrt(k[i] * k[i] + k[m] * k[m] - (2. * k[i] * k[m] * gamma_nm));

			double Rm = k[m] * tanh(k[m] * depth);

			double Rsm = sqrt(Rm);
			
			double Rsnm2plus = pow(Rsn + Rsm, 2.);
			double Rsnm2minus = pow(Rsn - Rsm, 2.);

			double D_nm_plus = (2. * Rsnm2plus * (k[i] * k[m] * gamma_nm - Rn * Rm)) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth)) + 
			((Rsn + Rsm) * (Rsn*(k[m]*k[m] - Rm * Rm) + Rsm*(k[i] * k[i] - Rn * Rn))) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth));

			double D_nm_minus = (2. * Rsnm2minus * (k[i] * k[m] * gamma_nm + Rn * Rm)) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth)) 
			+ ((Rsn - Rsm) * (-Rsn*(k[m]*k[m] - Rm * Rm) + Rsm*(k[i]*k[i] - Rn * Rn))) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth));

			double AA = 0.5*(A[i] * g / omega[i]) * (A[m] * g / omega[m]);
			double BB = (D_nm_plus/(omega[i]+omega[m]));
			double coshcosh_plus = cosh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth);
			
			// Sum terms
			vsum2 += AA * BB * (k[i]*sin(thi) + k[m]*sin(thm)) * cos(phi_i + phi_m) * coshcosh_plus;
	
			if (k[i] != k[m]){
				double CC = (D_nm_minus/(omega[i]-omega[m])); 
				double coshcosh_minus = cosh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth);
				
				// Diff terms
				vsum2 += AA * CC * (k[i]*sin(thi)-k[m]*sin(thm)) * cos(phi_i-phi_m) * coshcosh_minus;
			}					
		}
		vsum2_t += vsum2*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);	
	}
	return vsum2_t;
}

/* Second order horizontal velocity V for a sinus wave */
double Irregular::v2_boost(double t, double xx, double yy, double zz) {

	//double eta1_t = 0;
	double vsum2_t = 0;
	#pragma omp parallel for reduction(+:vsum2_t)
	for (int i = 0; i < nfreq * ndir; i++) {
		double vsum2 = 0;		
		double thi = theta[i]+ (mtheta * PI / 180.);
		double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];
	
		// Diagonal terms
		vsum2 += uvw2_diag[i] * sin(thi) * cos(2*phi_i) * cosh(2*k[i]*(zz+depth)) / cosh(2*k[i]*depth); 

		for (int m = i + 1; m < nfreq * ndir; m++) {
			double thm = theta[m]+ (mtheta * PI / 180.);
			
			double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];

			double gamma_nm = cos(thi - thm);
			double k_nm_plus = sqrt(k[i] * k[i] + k[m] * k[m] + (2. * k[i] * k[m] * gamma_nm));
			double k_nm_minus = sqrt(k[i] * k[i] + k[m] * k[m] - (2. * k[i] * k[m] * gamma_nm));

			double coshcosh_plus = cosh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth);
			double coshcosh_minus = cosh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth);
			
			vsum2 += uvw2_sum[i*(nfreq * ndir)+ m] * (k[i]*sin(thi) + k[m]*sin(thm)) * cos(phi_i + phi_m) * coshcosh_plus;
			vsum2 += uvw2_diff[i*(nfreq * ndir)+ m] * (k[i]*sin(thi) - k[m]*sin(thm)) * cos(phi_i - phi_m) * coshcosh_minus;					
		}
		vsum2_t += vsum2*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);	
	}
	return vsum2_t;
}

/* Second order vertical velocity component W for a sinus wave */
double Irregular::w2(double t, double xx, double yy, double zz) {

	//double eta1_t = 0;
	double wsum2_t = 0;
	//double kmax = *max_element(k.begin(), k.end());
	//double kmin = *min_element(k.begin(), k.end());
	for (int i = 0; i < nfreq * ndir; i++) {
		double wsum2 = 0;
		// Second order
		double thi = theta[i]+ (mtheta * PI / 180.);
		double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];
		double Rn = k[i] * tanh(k[i] * depth);
		double Rsn = sqrt(Rn);

		// Diagonal terms
		double D_nn_plus = ( 6.*Rsn*Rsn * (k[i] * k[i] - Rn * Rn))/(2.* Rsn* Rsn - k[i] * tanh(2*k[i] * depth));		
		double AA = 0.25*(A[i] * g/ omega[i])*(A[i] * g / omega[i])*(D_nn_plus/(2*omega[i]));
		wsum2 += AA* 2*k[i] * sin(2*phi_i) * sinh(2*k[i]*(zz+depth)) / cosh(2*k[i]*depth); 

		for (int m = i + 1; m < nfreq * ndir; m++) {
			double thm = theta[m]+ (mtheta * PI / 180.);
			
			double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];

			double gamma_nm = cos(thi - thm);
			double k_nm_plus = sqrt(k[i] * k[i] + k[m] * k[m] + (2. * k[i] * k[m] * gamma_nm));
			double k_nm_minus = sqrt(k[i] * k[i] + k[m] * k[m] - (2. * k[i] * k[m] * gamma_nm));

			double Rm = k[m] * tanh(k[m] * depth);

			double Rsm = sqrt(Rm);
			
			double Rsnm2plus = pow(Rsn + Rsm, 2.);
			double Rsnm2minus = pow(Rsn - Rsm, 2.);

			double D_nm_plus = (2. * Rsnm2plus * (k[i] * k[m] * gamma_nm - Rn * Rm)) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth)) + 
			((Rsn + Rsm) * (Rsn*(k[m]*k[m] - Rm * Rm) + Rsm*(k[i] * k[i] - Rn * Rn))) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth));

			double D_nm_minus = (2. * Rsnm2minus * (k[i] * k[m] * gamma_nm + Rn * Rm)) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth)) + ((Rsn - Rsm) * (-Rsn*(k[m]*k[m] - Rm * Rm) + Rsm*(k[i]*k[i] - Rn * Rn))) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth));

			double AA = 0.5*(A[i] * g / omega[i]) * (A[m] * g / omega[m]);
			double BB = (D_nm_plus/(omega[i]+omega[m]));
			
			double sinhcosh_plus = sinh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth);

			// Sum term
			wsum2 += AA * BB * k_nm_plus * sin(phi_i+phi_m) * sinhcosh_plus;	
			if (k[i] != k[m]){
				double CC = (D_nm_minus/(omega[i]-omega[m])); 
				double sinhcosh_minus = sinh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth);				
				// Diff terms
				wsum2 += AA * CC * k_nm_minus * sin(phi_i-phi_m) * sinhcosh_minus;
			}
		}
		wsum2_t += wsum2*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return wsum2_t;
}

/* Second order vertical velocity component W for a sinus wave */
double Irregular::w2_boost(double t, double xx, double yy, double zz) {

	//double eta1_t = 0;
	double wsum2_t = 0;
	#pragma omp parallel for reduction(+:wsum2_t)
	for (int i = 0; i < nfreq * ndir; i++) {
		double wsum2 = 0;
		// Second order
		double thi = theta[i]+ (mtheta * PI / 180.);
		double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];
		

		// Diagonal terms
		
		wsum2 += uvw2_diag[i] * sin(2*phi_i) * sinh(2*k[i]*(zz+depth)) / cosh(2*k[i]*depth); 

		for (int m = i + 1; m < nfreq * ndir; m++) {
			double thm = theta[m]+ (mtheta * PI / 180.);
			double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];

			double gamma_nm = cos(thi - thm);
			double k_nm_plus = sqrt(k[i] * k[i] + k[m] * k[m] + (2. * k[i] * k[m] * gamma_nm));
			double k_nm_minus = sqrt(k[i] * k[i] + k[m] * k[m] - (2. * k[i] * k[m] * gamma_nm));

			double sinhcosh_plus = sinh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth);
			double sinhcosh_minus = sinh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth);	

			// Sum term
			wsum2 += uvw2_sum[i*(nfreq * ndir)+ m] * k_nm_plus * sin(phi_i+phi_m) * sinhcosh_plus;	
				
			// Diff terms
			wsum2 += uvw2_diff[i*(nfreq * ndir)+ m] * k_nm_minus * sin(phi_i-phi_m) * sinhcosh_minus;
			
		}
		wsum2_t += wsum2*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return wsum2_t;
}


/* Second order dynamic pressure */
double Irregular::dp2(double t, double xx, double yy, double zz) {

	double psum2_t = 0;
	
	for (int i = 0; i < nfreq * ndir; i++) {
		double psum2 = 0;
	
		double thi = theta[i]+ (mtheta * PI / 180.);
		double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];
		double Rn = k[i] * tanh(k[i] * depth);
		double Rsn = sqrt(Rn);

		// Diagonal terms
		double D_nn_plus = ( 6.*Rsn*Rsn * (k[i] * k[i] - Rn * Rn))/(2.* Rsn* Rsn - k[i] * tanh(2*k[i] * depth));		
		double AA = 0.25*(A[i] * g/ omega[i])*(A[i] * g / omega[i])*(D_nn_plus/(2*omega[i]));
		psum2 += AA*(2*omega[i]) * cos(2*phi_i) * cosh(2*k[i]*(zz+depth)) / cosh(2*k[i]*depth); // first term
		psum2 -= (g*g*(2*k[i]*cos(thi)+2*k[i]*sin(thi))*2*profileX(i, xx, yy, zz))/(2*omega[i]); // second term
		psum2 += (g*g*(2*k[i])*2*profileZ(i, xx, yy, zz))/(2*omega[i]); // second term

		// Off-diagonals
		for (int m = i + 1; m < nfreq * ndir; m++) {
			double thm = theta[m]+ (mtheta * PI / 180.);
			
			double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];

			double gamma_nm = cos(thi - thm);
			double k_nm_plus = sqrt(k[i] * k[i] + k[m] * k[m] + (2. * k[i] * k[m] * gamma_nm));
			double k_nm_minus = sqrt(k[i] * k[i] + k[m] * k[m] - (2. * k[i] * k[m] * gamma_nm));

			double Rm = k[m] * tanh(k[m] * depth);

			double Rsm = sqrt(Rm);
			
			double Rsnm2plus = pow(Rsn + Rsm, 2.);
			double Rsnm2minus = pow(Rsn - Rsm, 2.);

			double D_nm_plus = (2. * Rsnm2plus * (k[i] * k[m] * gamma_nm - Rn * Rm)) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth)) + 
			((Rsn + Rsm) * (Rsn*(k[m]*k[m] - Rm * Rm) + Rsm*(k[i] * k[i] - Rn * Rn))) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth));

			double D_nm_minus = (2. * Rsnm2minus * (k[i] * k[m] * gamma_nm + Rn * Rm)) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth)) + ((Rsn - Rsm) * (-Rsn*(k[m]*k[m] - Rm * Rm) + Rsm*(k[i]*k[i] - Rn * Rn))) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth));

			double AA = 0.5*(A[i] * g / omega[i]) * (A[m] * g / omega[m]);
			double BB = (D_nm_plus/(omega[i]+omega[m]));
			double coshcosh_plus = cosh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth);

			// Sum terms
			psum2 += AA * BB * (omega[i] + omega[m]) * cos(phi_i + phi_m) * coshcosh_plus; //first term
			psum2 -= (g*g*(k[i]*cos(thi)*k[m]*cos(thm)+k[i]*sin(thi)*k[m]*sin(thm))*profileX(i, xx, yy, zz)*profileX(m, xx, yy, zz))/(omega[i]*omega[m]); // second term
			psum2 += (g*g*(k[i]*k[m])*profileZ(i, xx, yy, zz)*profileZ(m, xx, yy, zz))/(omega[i]*omega[m]); // second term
	
			if (k[i] != k[m]){
				double CC = (D_nm_minus/(omega[i]-omega[m])); 
				double coshcosh_minus = cosh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth);
				
				// Diff terms
				psum2 += AA * CC * (omega[i]-omega[m]) * cos(phi_i-phi_m) * coshcosh_minus;
			}	
		}
		psum2_t += psum2*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	// dp2 terms above were assembled without the density factor; dp1 includes rho,
	// so multiply here to keep the second-order dynamic pressure dimensionally
	// consistent (dynamic pressure = rho * [potential + kinetic terms]).
	return rho * psum2_t;
}

void Irregular::rebuild_component_cache()
{
    const int N = nfreq * ndir;
    const double mt = mtheta * PI / 180.0;

    cth_cache.resize(N);
    sth_cache.resize(N);
    R_cache.resize(N);
    Rs_cache.resize(N);
    ag_over_om_cache.resize(N);
    cosh2kd_cache.resize(N);
    phase0_cache.resize(N);

    kmax_cache = 0.0;

    for (int i = 0; i < N; ++i) {
        const double th = theta[i] + mt;

        cth_cache[i] = std::cos(th);
        sth_cache[i] = std::sin(th);

        const double R = k[i] * std::tanh(k[i] * depth);
        R_cache[i] = R;
        Rs_cache[i] = std::sqrt(R);

        ag_over_om_cache[i] = A[i] * g / omega[i];
        cosh2kd_cache[i] = std::cosh(2.0 * k[i] * depth);
        phase0_cache[i] = omega[i] * tofmax + phase[i];

        if (k[i] > kmax_cache) {
            kmax_cache = k[i];
        }
    }
}

void Irregular::rebuild_pair_cache()
{
    const int N = nfreq * ndir;
    pair_cache.clear();
    pair_cache.reserve(N * (N - 1) / 2);

    for (int i = 0; i < N; ++i) {
        if (omega[i] > dw_cutoff) {
            continue;
        }

        const double ki = k[i];
        const double oi = omega[i];
        const double Rni = R_cache[i];
        const double Rsni = Rs_cache[i];

        for (int m = i + 1; m < N; ++m) {
            if (omega[m] > dw_cutoff) {
                continue;
            }
            if (std::abs(oi - omega[m]) > dw_bandwidth) {
                continue;
            }

            const double km = k[m];
            const double om = omega[m];
            const double Rnm = R_cache[m];
            const double Rsm = Rs_cache[m];

            const double gamma_nm =
                cth_cache[i] * cth_cache[m] +
                sth_cache[i] * sth_cache[m];

            const double ki2 = ki * ki;
            const double km2 = km * km;
            const double kikm = ki * km;

            const double k_nm_plus =
                std::sqrt(ki2 + km2 + 2.0 * kikm * gamma_nm);

            const double k_nm_minus =
                std::sqrt(ki2 + km2 - 2.0 * kikm * gamma_nm);

            const double sp = Rsni + Rsm;
            const double sm = Rsni - Rsm;
            const double Rsnm2plus = sp * sp;
            const double Rsnm2minus = sm * sm;

            PairCacheEntry entry{};
            entry.i = i;
            entry.m = m;
            entry.gamma_nm = gamma_nm;
            entry.k_nm_plus = k_nm_plus;
            entry.k_nm_minus = k_nm_minus;
            entry.has_plus = false;
            entry.has_minus = false;
            entry.D_nm_plus = 0.0;
            entry.D_nm_minus = 0.0;
            entry.alpha_nm_plus = 0.0;
            entry.alpha_nm_minus = 0.0;
            entry.H_nm_plus = 0.0;
            entry.H_nm_minus = 0.0;

            if ((ki + km) <= kmax_cache) {
                const double denom_plus =
                    Rsnm2plus - k_nm_plus * std::tanh(k_nm_plus * depth);

                const double D_nm_plus =
                    (2.0 * Rsnm2plus * (kikm * gamma_nm - Rni * Rnm)) / denom_plus +
                    (sp * (Rsni * (km2 - Rnm * Rnm) + Rsm * (ki2 - Rni * Rni))) / denom_plus;

                const double alpha_nm_plus =
                    (D_nm_plus - (kikm * gamma_nm - Rni * Rnm)) / std::sqrt(Rni * Rnm)
                    + (Rni + Rnm);

                entry.has_plus = true;
                entry.D_nm_plus = D_nm_plus;
                entry.alpha_nm_plus = alpha_nm_plus;

                // Smit et al. (2017) H coefficient for sum interaction
                // Derived by inverting Eq. 10: C = (iω₁₂/g)H + ... with C = α/2
                // H = i·H_real, where H_real = -g/ω₁₂ · [α/2 - (ω₁²+ω₂²+ω₁ω₂)/(2g) + g·k₁·k₂/(2ω₁ω₂)]
                {
                    const double om12 = oi + om;
                    const double k1dk2 = kikm * gamma_nm;
                    entry.H_nm_plus = -g / om12 * (
                        alpha_nm_plus / 2.0
                        - (oi * oi + om * om + oi * om) / (2.0 * g)
                        + g * k1dk2 / (2.0 * oi * om));
                }
            }

            if (ki != km) {
                const double denom_minus =
                    Rsnm2minus - k_nm_minus * std::tanh(k_nm_minus * depth);

                const double D_nm_minus =
                    (2.0 * Rsnm2minus * (kikm * gamma_nm + Rni * Rnm)) / denom_minus +
                    (sm * (-Rsni * (km2 - Rnm * Rnm) + Rsm * (ki2 - Rni * Rni))) / denom_minus;

                const double alpha_nm_minus =
                    (D_nm_minus - (kikm * gamma_nm + Rni * Rnm)) / std::sqrt(Rni * Rnm)
                    + (Rni + Rnm);

                entry.has_minus = true;
                entry.D_nm_minus = D_nm_minus;
                entry.alpha_nm_minus = alpha_nm_minus;

                // Smit H coefficient for difference interaction
                // Derived by inverting Eq. 10 with ω₂→−ω₂, k₂→−k₂
                {
                    const double om12 = oi - om;
                    const double k1dk2 = kikm * gamma_nm;
                    if (std::abs(om12) > SMALL) {
                        entry.H_nm_minus = -g / om12 * (
                            alpha_nm_minus / 2.0
                            - (oi * oi + om * om - oi * om) / (2.0 * g)
                            + g * k1dk2 / (2.0 * oi * om));
                    }
                }
            }

            if (entry.has_plus || entry.has_minus) {
                pair_cache.push_back(entry);
            }
        }
    }
}

void Irregular::rebuild_second_order_cache()
{
    rebuild_component_cache();
    rebuild_pair_cache();
}

UVWP Irregular::uvw2(double t, double xx, double yy, double zz)
{
    const int N = nfreq * ndir;
    const double dxp = xx - fpoint[0];
    const double dyp = yy - fpoint[1];
    const double zpd = zz + depth;
    const double gg = g * g;

    if ((int)cth_cache.size() != N) {
        rebuild_second_order_cache();
    }

    double usum2_t = 0.0;
    double vsum2_t = 0.0;
    double wsum2_t = 0.0;
    double psum2_t = 0.0;

    std::vector<double> cphi(N), sphi(N), ramp(N), px(N), pz(N);

    for (int i = 0; i < N; ++i) {
        const double ph =
            k[i] * (cth_cache[i] * dxp + sth_cache[i] * dyp)
            - omega[i] * t
            + phase0_cache[i];

        cphi[i] = std::cos(ph);
        sphi[i] = std::sin(ph);
        ramp[i] = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
        px[i] = profileX(i, xx, yy, zz);
        pz[i] = profileZ(i, xx, yy, zz);
    }

    // diagonal terms
    for (int i = 0; i < N; ++i) {
        if (omega[i] > dw_cutoff) {
            continue;
        }

        const double ki = k[i];
        const double oi = omega[i];
        const double cti = cth_cache[i];
        const double sti = sth_cache[i];
        const double Rni = R_cache[i];
        const double Rsni = Rs_cache[i];
        const double agi = ag_over_om_cache[i];

        const double D_nn_plus =
            (6.0 * Rsni * Rsni * (ki * ki - Rni * Rni)) /
            (2.0 * Rsni * Rsni - ki * std::tanh(2.0 * ki * depth));

        const double AA_diag = 0.25 * agi * agi * (D_nn_plus / (2.0 * oi));

        const double cos_2phi_i = cphi[i] * cphi[i] - sphi[i] * sphi[i];
        const double sin_2phi_i = 2.0 * sphi[i] * cphi[i];

        const double coshcosh_diag = std::cosh(2.0 * ki * zpd) / cosh2kd_cache[i];
        const double sinhcosh_diag = std::sinh(2.0 * ki * zpd) / cosh2kd_cache[i];

        double u = 0.0, v = 0.0, w = 0.0, p = 0.0;

        u += AA_diag * (2.0 * ki * cti) * cos_2phi_i * coshcosh_diag;
        v += AA_diag * (2.0 * ki * sti) * cos_2phi_i * coshcosh_diag;
        w += AA_diag * (2.0 * ki) * sin_2phi_i * sinhcosh_diag;
        p += AA_diag * (2.0 * oi) * cos_2phi_i * coshcosh_diag;

        p -= (gg * (2.0 * ki * cti + 2.0 * ki * sti) * 2.0 * px[i]) / (2.0 * oi);
        p += (gg * (2.0 * ki) * 2.0 * pz[i]) / (2.0 * oi);

        const double ri = ramp[i];
        usum2_t += u * ri;
        vsum2_t += v * ri;
        wsum2_t += w * ri;
        psum2_t += p * ri;
    }

    // pair terms
    for (const auto& pair : pair_cache) {
        const int i = pair.i;
        const int m = pair.m;

        const double ki = k[i];
        const double km = k[m];
        const double oi = omega[i];
        const double om = omega[m];
        const double cti = cth_cache[i];
        const double sti = sth_cache[i];
        const double ctm = cth_cache[m];
        const double stm = sth_cache[m];
        const double ri = ramp[i];

        const double AA = 0.5 * ag_over_om_cache[i] * ag_over_om_cache[m];

        if (pair.has_plus) {
            const double BB = pair.D_nm_plus / (oi + om);
            const double cos_phi_plus = cphi[i] * cphi[m] - sphi[i] * sphi[m];
            const double sin_phi_plus = sphi[i] * cphi[m] + cphi[i] * sphi[m];

            const double denom = std::cosh(pair.k_nm_plus * depth);
            const double coshcosh_plus = std::cosh(pair.k_nm_plus * zpd) / denom;
            const double sinhcosh_plus = std::sinh(pair.k_nm_plus * zpd) / denom;

            usum2_t += AA * BB * (ki * cti + km * ctm) * cos_phi_plus * coshcosh_plus * ri;
            vsum2_t += AA * BB * (ki * sti + km * stm) * cos_phi_plus * coshcosh_plus * ri;
            wsum2_t += AA * BB * pair.k_nm_plus * sin_phi_plus * sinhcosh_plus * ri;
            psum2_t += AA * BB * (oi + om) * cos_phi_plus * coshcosh_plus * ri;

            psum2_t -= (gg * (ki * cti * km * ctm + ki * sti * km * stm) * px[i] * px[m]) / (oi * om) * ri;
            psum2_t += (gg * (ki * km) * pz[i] * pz[m]) / (oi * om) * ri;
        }

        if (pair.has_minus) {
            const double CC = pair.D_nm_minus / (oi - om);
            const double cos_phi_minus = cphi[i] * cphi[m] + sphi[i] * sphi[m];
            const double sin_phi_minus = sphi[i] * cphi[m] - cphi[i] * sphi[m];

            const double denom = std::cosh(pair.k_nm_minus * depth);
            const double coshcosh_minus = std::cosh(pair.k_nm_minus * zpd) / denom;
            const double sinhcosh_minus = std::sinh(pair.k_nm_minus * zpd) / denom;

            usum2_t += AA * CC * (ki * cti - km * ctm) * cos_phi_minus * coshcosh_minus * ri;
            vsum2_t += AA * CC * (ki * sti - km * stm) * cos_phi_minus * coshcosh_minus * ri;
            wsum2_t += AA * CC * pair.k_nm_minus * sin_phi_minus * sinhcosh_minus * ri;
            psum2_t += AA * CC * (oi - om) * cos_phi_minus * coshcosh_minus * ri;
        }
    }

    return {usum2_t, vsum2_t, wsum2_t, psum2_t};
}




double Irregular::dp1(double t, double xx, double yy, double zz) {

	double psum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		psum += A[i] * rho * g * (cosh(k[i] * (zz + depth)) / cosh(k[i] * depth)) 
			* cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[i] 
				+ (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi)*
				freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}

	return psum;

}

/* vertical velocity gradient at z=0 for velocity component U */
double Irregular::phi1_dxdz(double t, double xx, double yy) {

	double usum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		usum += cos(theta[i] + (mtheta * PI / 180.)) * A[i] * omega[i] * k[i] 
			* cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[i] 
				+ (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi)*
				freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return usum;

}

/* vertical velocity gradient at z=0 for velocity component U */
double Irregular::dp1_dz(double t, double xx, double yy) {

	double psum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		psum +=  A[i] * rho * g * k[i] 
			* cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[i] 
				+ (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi)*
				freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return psum;

}


/* vertical velocity gradient at z=0 for velocity component V */
double Irregular::phi1_dydz(double t, double xx, double yy) {

	double vsum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		vsum += sin(theta[i] + (mtheta * PI / 180.)) * A[i] * omega[i] * k[i] 
			* cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[i] 
				+ (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi)*
				freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return vsum;

}


/* vertical velocity gradient at z=0 for velocity component W */
double Irregular::phi1_dzdz(double t, double xx, double yy) {

	double wsum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		wsum += A[i] * omega[i] * k[i] * (cosh(k[i] * depth) / sinh(k[i] 
			* depth)) * sin(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) 
				+ sin(theta[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi)*
				freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return wsum;

}

// First order velocity potential
double Irregular::phi1_pot(double t, double xx, double yy, double zz) {

	double phisum = 0.0;
	double phi;


	for (int i = 0; i< ndir*nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		phisum += A[i] * (omega[i] / k[i]) * (cosh(k[i] * (zz + depth)) / sinh(k[i] * depth))*sin(k[i] 
		* (cos(theta[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(theta[i] + (mtheta*PI / 180.))*
		(yy - fpoint[1])) - omega[i] * t + phi)*
				freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}

	return phisum;

}


/* Second order component of velocity potential
Warning: Unsure if this is finished...should be checked thoroughly at some point*/
double Irregular::phi2_pot(double t, double xx, double yy, double zz) {
	//double eta1_t = 0;
	double phisum2_t = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		double phisum2 = 0;
		for (int j = 0; j < ndir; j++) {
			// Second order
			int ci = i * ndir + j;
			double phi_i = k[ci] * (cos(theta[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[ci] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			//cout << "Rn: " << Rn << endl;

			//// Adiusting Bandwidth for 2 order cut-off
			//if (i + 1 + f_bw<nfreqs) {
			//	p = i + 1 + f_bw;
			//}
			//else { p = nfreqs; }
			for (int m = i + 1; m < nfreq; m++) {
				int cm = m * ndir + j;
				if (abs(omega[cm] - omega[ci]) <= dw_bandwidth) {
					double gamma_nm = cos(theta[ci] - theta[cm]);
					double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2. * k[ci] * k[cm] * gamma_nm));
					double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2. * k[ci] * k[cm] * gamma_nm));

					double Rm = k[cm] * tanh(k[cm] * depth);


					double D_nm_plus = (sqrt(Rn) + sqrt(Rm)) * (sqrt(Rm) * (k[ci] * k[ci] - Rn * Rn) + sqrt(Rn) * (k[cm] * k[cm] - Rm * Rm)) + 2. * pow(sqrt(Rn) + sqrt(Rm), 2.) * (k[ci] * k[cm] * gamma_nm - Rn * Rm) /
						(pow(sqrt(Rn) + sqrt(Rm), 2.) - k_nm_plus * tanh(k_nm_plus * depth));
					double D_nm_minus = (sqrt(Rn) - sqrt(Rm)) * (sqrt(Rm) * (k[ci] * k[ci] - Rn * Rn) - sqrt(Rn) * (k[cm] * k[cm] - Rm * Rm)) + 2. * pow(sqrt(Rn) - sqrt(Rm), 2.) * (k[ci] * k[cm] * gamma_nm + Rn * Rm) /
						(pow(sqrt(Rn) - sqrt(Rm), 2.) - k_nm_minus * tanh(k_nm_minus * depth));


					double beta_nm_minus = D_nm_minus / (2 * k[ci] * k[cm] * (omega[ci] - omega[cm]));
					double beta_nm_plus = D_nm_plus / (2 * k[ci] * k[cm] * (omega[ci] + omega[cm]));

					// Catch NaN when two equal frequency components interact
					if (omega[ci] == omega[cm]) {
						D_nm_minus = 0.;
						beta_nm_minus = 0.;
					}

					double phi_m = k[cm] * (cos(theta[cm] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[cm] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[cm] * t + omega[cm] * tofmax + phase[cm];

					phisum2 += (A[ci] * A[cm] * omega[ci] * omega[cm]) *
						(beta_nm_minus * cos(phi_i - phi_m) * (k[ci] * cos(theta[ci] + (mtheta * PI / 180.)) - k[cm] * cos(theta[cm] + (mtheta * PI / 180.))) * (cosh(k_nm_minus * (zz + depth)) / cosh(k_nm_minus * depth)) +
							beta_nm_plus * cos(phi_i + phi_m) * (k[ci] * cos(theta[ci] + (mtheta * PI / 180.)) + k[cm] * cos(theta[cm] + (mtheta * PI / 180.))) * (cosh(k_nm_plus * (zz + depth)) / cosh(k_nm_plus * depth)));


				}
			}
		}
		phisum2_t += phisum2*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return phisum2_t;
}

double Irregular::sum(double ll[], int nsum) {
	double ss = 0.0;
	for (int i = 0; i < nsum; i++) {
		ss += ll[i];
	}
	return ss;
}

void Irregular::normalize_data() {
	// Normalize and/or amplify the amplitude spectrum if Normalize is switched on
	double sumA = std::accumulate(A.begin(), A.end(), 0.);
	if (normalize) {
		for (int i = 0; i < nfreq*ndir; i++) {
			A[i] = ampl * A[i] / sumA;
			phase[i] = 0.;
		}
	}
	else {
		for (int i = 0; i < nfreq*ndir; i++) {
			A[i] = ampl * A[i];
		}
	}
}


// trapezoidal integration function
double Irregular::trapz(double x[], double y[], int n) {
	
	double ss = 0.;

	for (int i = 0; i < (n-1); i++) {
		ss += (x[i+1]-x[i])*((y[i+1]+y[i])*0.5);
	}

	return ss;
}

double Irregular::interpolate(std::vector<double>& xData, std::vector<double>& yData, double x, bool extrapolate)
{
	int size = xData.size();

	int i = 0;                                                                  // find left end of interval for interpolation
	if (x >= xData[size - 2])                                                 // special case: beyond right end
	{
		i = size - 2;
	}
	else
	{
		while (x > xData[i + 1]) i++;
	}
	double xL = xData[i], yL = yData[i], xR = xData[i + 1], yR = yData[i + 1];      // points on either side (unless beyond ends)
	if (!extrapolate)                                                         // if beyond ends of array and not extrapolating
	{
		if (x < xL) yR = yL;
		if (x > xR) yL = yR;
	}

	double dydx = (yR - yL) / (xR - xL);                                    // gradient

	return yL + dydx * (x - xL);                                              // linear interpolation
}

std::vector<double>  Irregular::spectral_moments(){
	// Computes the spectral moments (m0, m1, m2) from spectral compnents
	if (ndir > 1) {
		std::cerr << "WARNING: mean wave length only supported for 1D wave spectra at the moment. Sorry for the inconvenience" << std::endl;
	}

	double w_t, dw, S;
	double m0 = 0.;
	double m1 = 0.;
	double m2 = 0.;

	// not very robust
	// for (int i = 0; i < (nfreq - 1); i++) {
	// 	w_t = (omega[i + 1] + omega[i])/2.;
	// 	dw = omega[i + 1] - omega[i];
	// 	S = pow((A[i + 1] + A[i]) * 0.5, 2.) / (2. * dw);
	// 	m0 += dw * S;
	// 	m1 += dw * S * w_t;
	// 	m2 += dw * S * w_t * w_t;
	// }

	for (int i = 0; i < (nfreq - 1); i++) {
		m0 += A[i]*A[i]/2.;
		m1 += omega[i]*A[i]*A[i]/2.;
		m2 += omega[i]*omega[i]*A[i]*A[i]/2.;
	}

	std::vector<double> res;
	res.push_back(m0);
	res.push_back(m1);
	res.push_back(m2);

	return res;
}

/* Function for calculating phase speed based integration of the wave spectrum*/
double Irregular::phase_speed(int opt)
{
	std::vector<double> moment = spectral_moments();
	
	double t1 = 2. * PI * (moment[0] / moment[1]);
	double t2 = 2. * PI * sqrt(moment[0] / moment[2]);

	double w1 = 2 * PI / t1;
	double w2 = 2 * PI / t2;

	// interpolate wave number
	double k1 = interpolate(omega, k, w1, true);
	double k2 = interpolate(omega, k, w2, true);

	if (opt == 1) {
		return w1 / k1;
	}
	else if(opt == 2) {
		return w2 / k2;
	}
	return 0.;
}




/* Function for calculating the mean wave length based integration of the wave spectrum*/
double Irregular::mean_wave_length(int opt)
{
	std::vector<double> moment = spectral_moments();
	double t1 = 2. * PI * (moment[0] / moment[1]);
	double t2 = 2. * PI * sqrt(moment[0] / moment[2]);

	double w1 = 2 * PI / t1;
	double w2 = 2 * PI / t2;

	// interpolate wave number
	double l1 = (2. * PI) / interpolate(omega, k, w1, true);
	double l2 = (2. * PI) / interpolate(omega, k, w2, true);

	if (opt == 1) {
		return l1;
	}
	else if (opt == 2) {
		return l2;
	}
	return 0.;
}

/* Function for calculating the mean wave length based integration of the wave spectrum*/
double Irregular::mean_wave_period(int opt)
{
	std::vector<double> moment = spectral_moments();

	double t1 = 2. * PI * (moment[0] / moment[1]);
	double t2 = 2. * PI * sqrt(moment[0] / moment[2]);

	if (opt == 1) {
		return t1;
	}
	else if (opt == 2) {
		return t2;
	}
	return 0.;
}

/* Function for calculating the mean wave length based integration of the wave spectrum*/
double Irregular::bandwidth_estimator()
{
	// check 
	if (ndir > 1) {
		std::cerr << "Automatic calculation of bandwidth parameter only supported for 2D wave spectra at the moment. Please specify bandiwidth manually. Sorry for the inconvenience" << std::endl;
		exit(-1);
	}

	double w_t, dw, S;
	std::vector<double> moment = spectral_moments();

	double bw = 0.7 * ( moment[1] / moment[0]);
	
	return bw;
}

// when called, writes spectral components stored in irregular class to dump file.
void Irregular::dumpSpectralComponents() {
	std::string fpath = "./spectral_components.dat";
	FILE* fp = fopen(fpath.c_str(), "w");
	fprintf(fp, "# spectral wave data output\n");
	fprintf(fp, "# [irregular wave components]\n");
	fprintf(fp, "# nfreq %d\n",nfreq);
	fprintf(fp, "# ndir 0\n");
	fprintf(fp, "# omega [rad/s]    a[m]       k[rad/m]        phi[rad]    theta[rad]\n");
	for (int i = 0; i < nfreq; i++) {
		fprintf(fp, "%12.5E  %12.5E  %12.5E  %12.5E  %12.5E\n", omega[i], A[i], k[i], phase[i], theta[i]);
	}
	fclose(fp);
}

void Irregular::print_spectral_parameters() {
	// Print some spectral parameters which can be used to validate the input spectral components
	// Spectral moments 0, 1, 2
	// Tm, Tz, Hs,
	if (A.empty() || omega.empty()) {
		std::cerr << "WARNING: no spectral components available -- skipping spectral parameter report. "
		             "Check that the [wave input irregular] section was parsed (correct header and frequency count)." << std::endl;
		return;
	}
	std::vector<double> moment = spectral_moments();
	double t1 = 2. * PI * (moment[0] / moment[1]);
	double t2 = 2. * PI * sqrt(moment[0] / moment[2]);
	double hs = 4*sqrt(moment[0]);
	std::vector<double>::iterator it = std::max_element(A.begin(), A.end());
    int index = static_cast<int>(std::distance(A.begin(), it));
	double tp = 2*PI/omega[index];
	

	std::cout << "----------------------------------------------------------------" << std::endl;
	std::cout << "Spectral parameters computed from specified spectral components:" << std::endl << std::endl;
	std::cout << "Spectral moments m0: " << moment[0] << ", m1: " << moment[1] << ", m2: " << moment[2] << std::endl;
	std::cout << "Significant Wave height Hs (m):" << hs << std::endl;
	std::cout << "Mean period Tz (sec, from m0 and m2):" << t2 << std::endl;
	std::cout << "Mean period Tm:(sec, from m0 and m1)" << t1 << std::endl;
	std::cout << "Spectral peak period (sec):" << tp << std::endl;
	std::cout << "Mean phase speed (m/s, based on Tm):" << phase_speed(1) << std::endl;
	std::cout << "----------------------------------------------------------------" << std::endl;



}


/* ==========================================================================
   Smit et al. (2017) s-coordinate second-order kinematics.

   The velocity field is (Eq. 11):
     u = Σ η̂₁ U₁¹ e₁  +  Σ η̂₁η̂₂ U¹² e₁e₂

   First-order uses Wheeler-stretched profiles: Ch^z = cosh[k·h·(h+z)/(h+η)] / cosh(kh)
   Second-order has the bound wave (H coefficient + stretched Ch_{1+2}) plus
   a coordinate-transformation correction term g(h+z)/(2(h+η)) · (...)
   ========================================================================== */

/* Second-order surface elevation built ENTIRELY from Smit's published formulas:
   H from Eq. 9 (paper formula, with the colleague's typo corrections applied),
   then C from Eq. 10, then eta2 = sum over pairs of C * eta_hat_i * eta_hat_m.

   This is a deliberately independent route to eta2 so it can be compared, term
   by term, against the trusted Sharma & Dean eta2(). If the two disagree we
   have isolated a transcription error in Smit's equations.

   Conventions (matching SD eta2()):
     - off-diagonal sum, per unordered pair {i,m}:  A_i A_m C_nm cos(phi_i+phi_m)
       with C = alpha/2, so this equals 0.5 A_i A_m alpha cos(...)  [== SD]
     - diagonal:  0.5 A_i^2 C_nn cos(2 phi_i)  =  0.25 A_i^2 alpha_nn cos(...) */
double Irregular::eta_smit(double t, double xx, double yy)
{
	const int N = nfreq * ndir;
	const double dxp = xx - fpoint[0];
	const double dyp = yy - fpoint[1];
	if ((int)cth_cache.size() != N) rebuild_second_order_cache();

	std::vector<double> cphi(N), sphi(N), ramp(N);
	for (int i = 0; i < N; ++i) {
		const double ph = k[i] * (cth_cache[i] * dxp + sth_cache[i] * dyp)
			- omega[i] * t + phase0_cache[i];
		cphi[i] = std::cos(ph);
		sphi[i] = std::sin(ph);
		ramp[i] = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
	}

	double eta2_total = 0.0;

	// --- Diagonal terms (component i with itself, sum frequency 2*omega_i) ---
	for (int i = 0; i < N; ++i) {
		if (omega[i] > dw_cutoff) continue;
		if (2.0 * k[i] > kmax_cache) continue;

		const double ki = k[i], oi = omega[i];
		const double k1dk2 = ki * ki;          // k_i . k_i (collinear self)
		const double om12 = 2.0 * oi;
		const double k12 = 2.0 * ki;           // |k_i + k_i|
		const double sigma2 = g * k12 * std::tanh(k12 * depth);

		// Smit Eq. 9 (typo-corrected), H = i * H_real:
		//   H_real = om12 g^2 / (sigma^2 - om12^2) *
		//            [ k1.k2/(o1 o2) - (o1 o2 + o1^2 + o2^2)/g^2
		//              + (k1^2 o2 + k2^2 o1)/(2 o1 o2 om12) ]
		const double bracket =
			k1dk2 / (oi * oi)
			- (oi * oi + oi * oi + oi * oi) / (g * g)
			+ (ki * ki * oi + ki * ki * oi) / (2.0 * oi * oi * om12);
		const double H_real = om12 * g * g / (sigma2 - om12 * om12) * bracket;

		// Eq. 10:  C = (i om12/g) H + (o1^2+o2^2+o1 o2)/(2g) - g k1.k2/(2 o1 o2)
		//   with H = i H_real  =>  (i om12/g) H = -om12 H_real / g
		const double C = -om12 * H_real / g
			+ (oi * oi + oi * oi + oi * oi) / (2.0 * g)
			- g * k1dk2 / (2.0 * oi * oi);

		const double cos_2phi = cphi[i] * cphi[i] - sphi[i] * sphi[i];
		eta2_total += 0.5 * A[i] * A[i] * C * cos_2phi * ramp[i];
	}

	// --- Off-diagonal terms (unordered pairs) ---
	for (const auto& p : pair_cache) {
		const int i = p.i, m = p.m;
		const double ki = k[i], km = k[m];
		const double oi = omega[i], om = omega[m];
		const double k1dk2 = ki * km * p.gamma_nm;     // k_i . k_m

		if (p.has_plus) {
			const double om12 = oi + om;
			const double k12 = p.k_nm_plus;
			const double sigma2 = g * k12 * std::tanh(k12 * depth);
			const double bracket =
				k1dk2 / (oi * om)
				- (oi * om + oi * oi + om * om) / (g * g)
				+ (ki * ki * om + km * km * oi) / (2.0 * oi * om * om12);
			const double H_real = om12 * g * g / (sigma2 - om12 * om12) * bracket;
			const double C = -om12 * H_real / g
				+ (oi * oi + om * om + oi * om) / (2.0 * g)
				- g * k1dk2 / (2.0 * oi * om);
			const double cos_phi_plus = cphi[i] * cphi[m] - sphi[i] * sphi[m];
			eta2_total += A[i] * A[m] * C * cos_phi_plus * ramp[i];
		}

		if (p.has_minus) {
			const double om12 = oi - om;
			const double k12 = p.k_nm_minus;
			const double sigma2 = g * k12 * std::tanh(k12 * depth);
			// difference: omega_2 -> -omega_2, k_2 -> -k_2 in Eq. 9 & Eq. 10
			if (std::abs(om12) > SMALL) {
				const double bracket =
					k1dk2 / (oi * om)                       // k1.(-k2)/(o1.(-o2)) = k1.k2/(o1 o2)
					+ (oi * om - oi * oi - om * om) / (g * g)
					+ (ki * ki * (-om) + km * km * oi) / (2.0 * oi * (-om) * om12);
				const double H_real = om12 * g * g / (sigma2 - om12 * om12) * bracket;
				const double C = -om12 * H_real / g
					+ (oi * oi + om * om - oi * om) / (2.0 * g)
					- g * k1dk2 / (2.0 * oi * om);
				const double cos_phi_minus = cphi[i] * cphi[m] + sphi[i] * sphi[m];
				eta2_total += A[i] * A[m] * C * cos_phi_minus * ramp[i];
			}
		}
	}

	return eta2_total;
}


double Irregular::u_smit(double t, double xx, double yy, double zz) {

	// Step 1: compute surface elevation for the stretch factor
	double eta_val = eta1(t, xx, yy);
	if (order >= 1) eta_val += eta2(t, xx, yy);

	// Step 2: clamp to water column
	if (zz > eta_val + swl) return 0.0;
	zz = std::max(-depth, zz - swl);

	// Precompute phases and first-order stretched profiles
	const int N = ndir * nfreq;
	if ((int)cth_cache.size() != N) rebuild_second_order_cache();

	// Step 3: first-order velocity with stretched profiles (Wheeler stretching)
	double u1_val = 0.0;
	for (int i = 0; i < N; i++) {
		double phi = k[i] * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- omega[i] * t + phase0_cache[i];
		double fr = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		u1_val += A[i] * g * k[i] / omega[i] * profileX_s(i, zz, eta_val)
			* cth_cache[i] * cos(phi) * fr;
	}

	// Step 4: second-order contributions (skip if order < 1, giving pure Wheeler stretching)
	double u2_val = 0.0;
	if (order >= 1) {

	// Accumulate ∂_s φ² so the Eulerian chain term can map the FULL potential
	// (φ¹+φ²). The φ² part is the necessary O(ε³) term: −(1+s)/D·η_x·∂_s φ².
	double dphi2ds = 0.0;

	// Diagonal terms (self-interaction: component i with itself)
	for (int i = 0; i < N; i++) {
		double phi_i = k[i] * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- omega[i] * t + phase0_cache[i];
		double fr = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		double ki = k[i], oi = omega[i];
		double AA = A[i] * A[i] * fr * fr;

		double k2 = 2.0 * ki; // combined wavenumber magnitude
		double om2 = 2.0 * oi; // combined frequency
		// Compute diagonal alpha (Sharma & Dean) then invert Eq. 10 for H
		double Rni = R_cache[i], Rsni = Rs_cache[i];
		double D_nn = (6.0*Rsni*Rsni*(ki*ki - Rni*Rni)) / (2.0*Rsni*Rsni - ki*std::tanh(2.0*ki*depth));
		double alpha_nn = (D_nn - (ki*ki - Rni*Rni)) / Rni + 2.0*Rni;
		double H_diag = -g / om2 * (alpha_nn / 2.0
			- (3.0 * oi * oi) / (2.0 * g)
			+ g * ki * ki / (2.0 * oi * oi));

		// Bound wave: -kx₁₂ · H · Ch(2k) · cos(2φ)  [vertically-Lagrangian: ∂_x* at constant s]
		double kx2 = 2.0 * ki * cth_cache[i];
		double Ch2 = profileX_s2(k2, zz, eta_val);
		u2_val += 0.5 * AA * H_diag * kx2 * (-cos(2.0 * phi_i)) * Ch2;

		// P-derivative ∂ₓP (particular-solution velocity) — coefficient validated
		// against exact Stokes at all depths (papers/sigma_to_euler_derivation.py).
		double sp1 = (depth + zz) / (depth + eta_val);              // (1+s)
		double g1s = g * sp1;                                       // g(1+s)
		double Sh_i = profileZ_s(i, zz, eta_val);
		u2_val += 0.5 * AA * kx2 * g1s * (ki * Sh_i / oi) * cos(2.0 * phi_i);

		// ∂_s φ² (diagonal): bound  φ=−0.5·AA·H·Ch2·sin2φ,  ∂_sCh2 = k2·h·Sh2
		//                    P-pot  φ= 0.5·AA·g(1+s)(k·Sh/ω)·sin2φ
		double Ch_i = profileX_s(i, zz, eta_val);
		double Sh2  = profileZ_s2(k2, zz, eta_val);
		dphi2ds += -0.5 * AA * H_diag * (k2 * depth * Sh2) * sin(2.0 * phi_i);
		dphi2ds +=  0.5 * AA * g * (ki * Sh_i / oi
			+ sp1 * ki * (ki * depth * Ch_i) / oi) * sin(2.0 * phi_i);
	}

	// Off-diagonal terms (pair interactions)
	for (const auto& p : pair_cache) {
		const int i = p.i;
		const int m = p.m;
		const double Ai = A[i], Am = A[m];
		const double ki = k[i], km = k[m];
		const double oi = omega[i], om = omega[m];

		double phi_i = ki * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- oi * t + phase0_cache[i];
		double phi_m = km * (cth_cache[m] * (xx - fpoint[0]) + sth_cache[m] * (yy - fpoint[1]))
			- om * t + phase0_cache[m];
		double fr_i = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		double fr_m = freq_ramp(t, freq_ramptime[m], freq_starttime[m], freq_duration[m]);
		double AA = Ai * Am * fr_i * fr_m;

		double sp1 = (depth + zz) / (depth + eta_val);                // (1+s)
		double g1s = g * sp1;                                         // g(1+s)
		double Sh_i = profileZ_s(i, zz, eta_val);
		double Sh_m = profileZ_s(m, zz, eta_val);
		double Ch_i = profileX_s(i, zz, eta_val);
		double Ch_m = profileX_s(m, zz, eta_val);

		// Sum interaction (ω₁ + ω₂): bound wave + P-derivative ∂ₓP
		if (p.has_plus) {
			double phi_plus = phi_i + phi_m;
			double kx12 = ki * cth_cache[i] + km * cth_cache[m];
			double Ch12 = profileX_s2(p.k_nm_plus, zz, eta_val);
			u2_val += AA * p.H_nm_plus * kx12 * (-cos(phi_plus)) * Ch12;
			u2_val += AA * kx12 * g1s * (ki * Sh_i / oi + km * Sh_m / om) * cos(phi_plus);

			// ∂_s φ² (sum): bound  φ=−AA·H·Ch12·sinφ₊,  P-pot φ=AA·g(1+s)(...)·sinφ₊
			double Sh12 = profileZ_s2(p.k_nm_plus, zz, eta_val);
			dphi2ds += -AA * p.H_nm_plus * (p.k_nm_plus * depth * Sh12) * sin(phi_plus);
			dphi2ds +=  AA * g * ((ki * Sh_i / oi + km * Sh_m / om)
				+ sp1 * (ki * (ki * depth * Ch_i) / oi + km * (km * depth * Ch_m) / om)) * sin(phi_plus);
		}

		// Difference interaction (ω₁ − ω₂): k₂→−k₂, ω₂→−ω₂
		if (p.has_minus) {
			double phi_minus = phi_i - phi_m;
			double kx12 = ki * cth_cache[i] - km * cth_cache[m];
			double Ch12 = profileX_s2(p.k_nm_minus, zz, eta_val);
			u2_val += AA * p.H_nm_minus * kx12 * (-cos(phi_minus)) * Ch12;
			u2_val += AA * kx12 * g1s * (ki * Sh_i / oi - km * Sh_m / om) * cos(phi_minus);

			// ∂_s φ² (difference): k₂→−k₂, ω₂→−ω₂
			double Sh12 = profileZ_s2(p.k_nm_minus, zz, eta_val);
			dphi2ds += -AA * p.H_nm_minus * (p.k_nm_minus * depth * Sh12) * sin(phi_minus);
			dphi2ds +=  AA * g * ((ki * Sh_i / oi - km * Sh_m / om)
				+ sp1 * (ki * (ki * depth * Ch_i) / oi - km * (km * depth * Ch_m) / om)) * sin(phi_minus);
		}
	}

	// Lagrangian→Eulerian chain term:  u += -(1+s)/D · η_x · φ_s
	// (field product of first-order surface slope η_x and ∂φ₁/∂s — automatically
	//  generates all sum/difference cross terms. Derived & verified symbolically in
	//  papers/sigma_to_euler_derivation.py: for one component this equals u1_quad,
	//  the missing half of the Eulerian 2nd-order velocity from φ₁.)
	{
		double eta_x = 0.0, phi_s = 0.0;
		for (int n = 0; n < N; n++) {
			double phi_n = k[n] * (cth_cache[n] * (xx - fpoint[0]) + sth_cache[n] * (yy - fpoint[1]))
				- omega[n] * t + phase0_cache[n];
			double fr = freq_ramp(t, freq_ramptime[n], freq_starttime[n], freq_duration[n]);
			double sphi = sin(phi_n);
			eta_x += -A[n] * k[n] * cth_cache[n] * sphi * fr;                       // ∂η₁/∂x
			phi_s += (g * A[n] / omega[n]) * k[n] * depth
				* profileZ_s(n, zz, eta_val) * sphi * fr;                           // ∂φ₁/∂s
		}
		phi_s += dphi2ds;   // include ∂_s φ²: chain maps the FULL potential (Eulerian, O(ε³) consistent)
		double oneplus_s_over_D = (depth + zz) / ((depth + eta_val) * (depth + eta_val));
		u2_val += -oneplus_s_over_D * eta_x * phi_s;
	}
	} // end if (order >= 1)

	return u1_val + u2_val;
}


double Irregular::v_smit(double t, double xx, double yy, double zz) {

	double eta_val = eta1(t, xx, yy);
	if (order >= 1) eta_val += eta2(t, xx, yy);
	if (zz > eta_val + swl) return 0.0;
	zz = std::max(-depth, zz - swl);

	const int N = ndir * nfreq;
	if ((int)cth_cache.size() != N) rebuild_second_order_cache();

	// First order with stretched profiles
	double v1_val = 0.0;
	for (int i = 0; i < N; i++) {
		double phi = k[i] * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- omega[i] * t + phase0_cache[i];
		double fr = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		v1_val += A[i] * g * k[i] / omega[i] * profileX_s(i, zz, eta_val)
			* sth_cache[i] * cos(phi) * fr;
	}

	// Second order
	double v2_val = 0.0;
	if (order >= 1) {

	// Diagonal terms
	for (int i = 0; i < N; i++) {
		double phi_i = k[i] * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- omega[i] * t + phase0_cache[i];
		double fr = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		double ki = k[i], oi = omega[i];
		double AA = A[i] * A[i] * fr * fr;
		double k2 = 2.0 * ki, om2 = 2.0 * oi;
		// Diagonal H from inverting Eq. 10 with Sharma & Dean alpha (NOT the buggy Eq. 9)
		double Rni = R_cache[i], Rsni = Rs_cache[i];
		double D_nn = (6.0*Rsni*Rsni*(ki*ki - Rni*Rni)) / (2.0*Rsni*Rsni - ki*std::tanh(2.0*ki*depth));
		double alpha_nn = (D_nn - (ki*ki - Rni*Rni)) / Rni + 2.0*Rni;
		double H_diag = -g / om2 * (alpha_nn / 2.0
			- (3.0 * oi * oi) / (2.0 * g)
			+ g * ki * ki / (2.0 * oi * oi));
		double ky2 = 2.0 * ki * sth_cache[i];
		double Ch2 = profileX_s2(k2, zz, eta_val);
		v2_val += 0.5 * AA * H_diag * ky2 * (-cos(2.0 * phi_i)) * Ch2;
	}

	// Off-diagonal terms
	for (const auto& p : pair_cache) {
		const int i = p.i, m = p.m;
		const double Ai = A[i], Am = A[m];
		const double ki = k[i], km = k[m];
		const double oi = omega[i], om = omega[m];

		double phi_i = ki * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- oi * t + phase0_cache[i];
		double phi_m = km * (cth_cache[m] * (xx - fpoint[0]) + sth_cache[m] * (yy - fpoint[1]))
			- om * t + phase0_cache[m];
		double fr_i = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		double fr_m = freq_ramp(t, freq_ramptime[m], freq_starttime[m], freq_duration[m]);
		double AA = Ai * Am * fr_i * fr_m;

		if (p.has_plus) {
			double phi_plus = phi_i + phi_m;
			double ky12 = ki * sth_cache[i] + km * sth_cache[m];
			double Ch12 = profileX_s2(p.k_nm_plus, zz, eta_val);
			v2_val += AA * p.H_nm_plus * ky12 * (-cos(phi_plus)) * Ch12;
		}

		if (p.has_minus) {
			double phi_minus = phi_i - phi_m;
			double ky12 = ki * sth_cache[i] - km * sth_cache[m];
			double Ch12 = profileX_s2(p.k_nm_minus, zz, eta_val);
			v2_val += AA * p.H_nm_minus * ky12 * (-cos(phi_minus)) * Ch12;
		}
	}
	} // end if (order >= 1)

	return v1_val + v2_val;
}


double Irregular::w_smit(double t, double xx, double yy, double zz) {

	double eta_val = eta1(t, xx, yy);
	if (order >= 1) eta_val += eta2(t, xx, yy);
	if (zz > eta_val + swl) return 0.0;
	zz = std::max(-depth, zz - swl);

	const int N = ndir * nfreq;
	if ((int)cth_cache.size() != N) rebuild_second_order_cache();

	// First order with stretched profiles
	double w1_val = 0.0;
	for (int i = 0; i < N; i++) {
		double phi = k[i] * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- omega[i] * t + phase0_cache[i];
		double fr = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		w1_val += A[i] * g * k[i] / omega[i] * profileZ_s(i, zz, eta_val)
			* sin(phi) * fr;
	}

	// Second order
	double w2_val = 0.0;
	if (order >= 1) {

	// Diagonal terms
	for (int i = 0; i < N; i++) {
		double phi_i = k[i] * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- omega[i] * t + phase0_cache[i];
		double fr = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		double ki = k[i], oi = omega[i];
		double AA = A[i] * A[i] * fr * fr;
		double k2 = 2.0 * ki, om2 = 2.0 * oi;
		// Diagonal H from inverting Eq. 10 with Sharma & Dean alpha (NOT the buggy Eq. 9)
		double Rni = R_cache[i], Rsni = Rs_cache[i];
		double D_nn = (6.0*Rsni*Rsni*(ki*ki - Rni*Rni)) / (2.0*Rsni*Rsni - ki*std::tanh(2.0*ki*depth));
		double alpha_nn = (D_nn - (ki*ki - Rni*Rni)) / Rni + 2.0*Rni;
		double H_diag = -g / om2 * (alpha_nn / 2.0
			- (3.0 * oi * oi) / (2.0 * g)
			+ g * ki * ki / (2.0 * oi * oi));
		// W bound: -|k₁₂|·H·Sh(2k)·sin(2φ)  (bound wave only, no correction term)
		double Sh2 = profileZ_s2(k2, zz, eta_val);
		w2_val += 0.5 * AA * H_diag * k2 * (-sin(2.0 * phi_i)) * Sh2;
	}

	// Off-diagonal terms
	for (const auto& p : pair_cache) {
		const int i = p.i, m = p.m;
		const double Ai = A[i], Am = A[m];
		const double ki = k[i], km = k[m];
		const double oi = omega[i], om = omega[m];

		double phi_i = ki * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- oi * t + phase0_cache[i];
		double phi_m = km * (cth_cache[m] * (xx - fpoint[0]) + sth_cache[m] * (yy - fpoint[1]))
			- om * t + phase0_cache[m];
		double fr_i = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		double fr_m = freq_ramp(t, freq_ramptime[m], freq_starttime[m], freq_duration[m]);
		double AA = Ai * Am * fr_i * fr_m;

		if (p.has_plus) {
			double phi_plus = phi_i + phi_m;
			double Sh12 = profileZ_s2(p.k_nm_plus, zz, eta_val);
			w2_val += AA * p.H_nm_plus * p.k_nm_plus * (-sin(phi_plus)) * Sh12;
		}

		if (p.has_minus) {
			double phi_minus = phi_i - phi_m;
			double Sh12 = profileZ_s2(p.k_nm_minus, zz, eta_val);
			w2_val += AA * p.H_nm_minus * p.k_nm_minus * (-sin(phi_minus)) * Sh12;
		}
	}
	} // end if (order >= 1)

	return w1_val + w2_val;
}


double Irregular::dp_smit(double t, double xx, double yy, double zz) {

	double eta_val = eta1(t, xx, yy);
	if (order >= 1) eta_val += eta2(t, xx, yy);
	if (zz > eta_val + swl) return 0.0;
	zz = std::max(-depth, zz - swl);

	const int N = ndir * nfreq;
	if ((int)cth_cache.size() != N) rebuild_second_order_cache();

	// First order pressure with stretched profiles
	double dp1_val = 0.0;
	for (int i = 0; i < N; i++) {
		double phi = k[i] * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- omega[i] * t + phase0_cache[i];
		double fr = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		dp1_val += A[i] * rho * g * profileX_s(i, zz, eta_val) * cos(phi) * fr;
	}

	// Second order pressure — bound wave + Bernoulli kinetic-energy term
	double dp2_val = 0.0;
	if (order >= 1) {

	// Diagonal pressure terms
	for (int i = 0; i < N; i++) {
		double phi_i = k[i] * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- omega[i] * t + phase0_cache[i];
		double fr = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		double ki = k[i], oi = omega[i];
		double AA = A[i] * A[i] * fr * fr;
		double k2 = 2.0 * ki, om2 = 2.0 * oi;
		// Diagonal H from inverting Eq. 10 with Sharma & Dean alpha (NOT the buggy Eq. 9)
		double Rni = R_cache[i], Rsni = Rs_cache[i];
		double D_nn = (6.0*Rsni*Rsni*(ki*ki - Rni*Rni)) / (2.0*Rsni*Rsni - ki*std::tanh(2.0*ki*depth));
		double alpha_nn = (D_nn - (ki*ki - Rni*Rni)) / Rni + 2.0*Rni;
		double H_diag = -g / om2 * (alpha_nn / 2.0
			- (3.0 * oi * oi) / (2.0 * g)
			+ g * ki * ki / (2.0 * oi * oi));
		// P bound: -ω₁₂·H·Ch(2k)·cos(2φ)·ρ
		double Ch2 = profileX_s2(k2, zz, eta_val);
		dp2_val += 0.5 * AA * rho * om2 * H_diag * (-cos(2.0 * phi_i)) * Ch2;
		// P kinetic (Bernoulli -rho/2 |grad phi|^2): -ρg²k²(Ch²−Sh²)/(ω²)·cos(2φ)
		double Ch_i = profileX_s(i, zz, eta_val);
		double Sh_i = profileZ_s(i, zz, eta_val);
		dp2_val -= 0.5 * AA * rho * g * g * ki * ki * (Ch_i * Ch_i - Sh_i * Sh_i) / (oi * oi) * cos(2.0 * phi_i);
	}

	// Off-diagonal pressure terms
	for (const auto& p : pair_cache) {
		const int i = p.i, m = p.m;
		const double Ai = A[i], Am = A[m];
		const double ki = k[i], km = k[m];
		const double oi = omega[i], om = omega[m];

		double phi_i = ki * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- oi * t + phase0_cache[i];
		double phi_m = km * (cth_cache[m] * (xx - fpoint[0]) + sth_cache[m] * (yy - fpoint[1]))
			- om * t + phase0_cache[m];
		double fr_i = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		double fr_m = freq_ramp(t, freq_ramptime[m], freq_starttime[m], freq_duration[m]);
		double AA = Ai * Am * fr_i * fr_m;

		if (p.has_plus) {
			double phi_plus = phi_i + phi_m;
			double om12 = oi + om;
			double Ch12 = profileX_s2(p.k_nm_plus, zz, eta_val);

			// Bound wave pressure: iω₁₂·H·Ch_{1+2} → ω₁₂·H·(−cos(phase))·ρ
			dp2_val += AA * rho * om12 * p.H_nm_plus * (-cos(phi_plus)) * Ch12;

			// Kinetic energy (Bernoulli): sum-freq bracket [Ch·Ch − Sh·Sh]
			double Ch_i = profileX_s(i, zz, eta_val);
			double Ch_m = profileX_s(m, zz, eta_val);
			double Sh_i = profileZ_s(i, zz, eta_val);
			double Sh_m = profileZ_s(m, zz, eta_val);
			double k1dk2 = ki * km * p.gamma_nm;
			dp2_val -= AA * rho * g * g * k1dk2 * (Ch_i * Ch_m - Sh_i * Sh_m) / (oi * om) * cos(phi_plus);
		}

		if (p.has_minus) {
			double phi_minus = phi_i - phi_m;
			double om12 = oi - om;
			double Ch12 = profileX_s2(p.k_nm_minus, zz, eta_val);

			dp2_val += AA * rho * om12 * p.H_nm_minus * (-cos(phi_minus)) * Ch12;

			// Kinetic energy (Bernoulli): diff-freq bracket [Ch·Ch + Sh·Sh]
			double Ch_i = profileX_s(i, zz, eta_val);
			double Ch_m = profileX_s(m, zz, eta_val);
			double Sh_i = profileZ_s(i, zz, eta_val);
			double Sh_m = profileZ_s(m, zz, eta_val);
			double k1dk2 = ki * km * p.gamma_nm;
			dp2_val -= AA * rho * g * g * k1dk2 * (Ch_i * Ch_m + Sh_i * Sh_m) / (oi * om) * cos(phi_minus);
		}
	}
	} // end if (order >= 1)

	return dp1_val + dp2_val;
}

