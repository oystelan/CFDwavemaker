#include "lsgrid_spline.h"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <string>
#include "Utils.h"
#if defined(_OPENMP)
#include "omp.h"
#else
// serial fallbacks for the omp timing/diagnostic calls, so the MPI build
// (which must NOT mix in OpenMP) compiles without libgomp
#include <chrono>
static inline double omp_get_wtime() {
	return std::chrono::duration<double>(
		std::chrono::steady_clock::now().time_since_epoch()).count();
}
static inline int omp_get_num_threads() { return 1; }
#endif
#include "mpi_utils.h"


//----------------------------------------------------------------------------------------------------------------------------------------


// LSgridSpline functions


//----------------------------------------------------------------------------------------------------------------------------------------



// Allocation of memory to storage matrices
void lsGridSpline::build_identity_slot() {
	// dense layout: storage index == i*ny + j, nothing dropped
	slot.resize((size_t)nx * ny);
	for (long k = 0; k < (long)nx * ny; k++) slot[k] = (int)k;
	col_off.resize(nx + 1);
	for (int i = 0; i <= nx; i++) col_off[i] = (long)i * ny;
	nkept = (long)nx * ny;
	trash = nkept;
	nnodes = nkept + 1;
	compacted = false;
}

/* Drop every node the ignore box covers, keeping a two-node halo around the
   retained set:

     ring 1  the bilinear / tricubic gather in square_vals()/cube_vals()
             reads the FAR corner (i+1, j+1) of the query cell, so a query
             anywhere in the retained area touches one node beyond it;
     ring 2  those corner nodes must carry VALID finite-difference
             gradients, and a gradient is built from the node's own
             neighbours;
     ring 3  margin, so that the gradient-valid region (kept minus one ring)
             extends TWO full nodes past the retained region rather than
             sitting exactly on its edge. Interpolation stays accurate for a
             client whose stencil reaches slightly beyond the box, instead of
             failing at the first node.

   Anything past ring 3 is never read, by either the interpolation or the
   gradient stencils, and is what this frees. One extra ring costs very
   little: the retained set is a thin frame, so its area grows roughly
   linearly with the halo width. Called either before the first
   fill (ignore_at_init: nothing to preserve) or after it (the retained data
   is carried across into the smaller arrays). */
void lsGridSpline::compact() {
	if (compacted) return;
	if (slot.size() != (size_t)nx * ny) build_identity_slot();

	const int HALO = 3;
	std::vector<char> keep((size_t)nx * ny, 0);
	long nignored = 0;
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++) {
			if (IGNORE[i * ny + j]) { nignored++; continue; }
			for (int a = -HALO; a <= HALO; a++) {
				const int ii = std::min(std::max(i + a, 0), nx - 1);
				for (int b = -HALO; b <= HALO; b++) {
					const int jj = std::min(std::max(j + b, 0), ny - 1);
					keep[(size_t)ii * ny + jj] = 1;
				}
			}
		}
	if (nignored == 0) return;   // no ignore box: nothing to gain

	// assign storage indices in i-major order, so each rank's column range
	// stays a CONTIGUOUS slot range (the MPI gather relies on that)
	const std::vector<int> old_slot = slot;
	const long old_nnodes = nnodes;
	std::vector<int> new_slot((size_t)nx * ny, -1);
	std::vector<long> new_col_off(nx + 1, 0);
	long k = 0;
	for (int i = 0; i < nx; i++) {
		new_col_off[i] = k;
		for (int j = 0; j < ny; j++)
			if (keep[(size_t)i * ny + j]) new_slot[(size_t)i * ny + j] = (int)k++;
	}
	new_col_off[nx] = k;
	const long new_nkept = k, new_nnodes = k + 1;

	if (new_nnodes >= old_nnodes) return;   // nothing to save

	// move each field into its smaller array, carrying the retained values
	struct Fld { double** p; int lev; bool layered; };
	std::vector<Fld> flds = {
		{&ETA,4,false},{&ETAdx,2,false},{&ETAdy,2,false},{&ETAdt,4,false},
		{&UX,4,true},{&UY,4,true},{&UZ,4,true},
		{&Udt,4,true},{&Udx,2,true},{&Udy,2,true},{&Uds,2,true},
		{&Vdt,4,true},{&Vdx,2,true},{&Vdy,2,true},{&Vds,2,true},
		{&Wdt,4,true},{&Wdx,2,true},{&Wdy,2,true},{&Wds,2,true}
	};
	if (store_pressure) {
		flds.push_back({&P,4,true});
		flds.push_back({&Pdt,4,true});
		flds.push_back({&Pdx,2,true});
		flds.push_back({&Pdy,2,true});
		flds.push_back({&Pds,2,true});
	}
	for (Fld& f : flds) {
		const long per = f.layered ? nl : 1;
		double* dst = new double[(size_t)f.lev * new_nnodes * per]();
		if (*f.p) {
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++) {
					const long ns = new_slot[(size_t)i * ny + j];
					const long os = old_slot[(size_t)i * ny + j];
					if (ns < 0 || os < 0) continue;
					for (int t = 0; t < f.lev; t++)
						for (long m = 0; m < per; m++)
							dst[(t * new_nnodes + ns) * per + m] =
								(*f.p)[(t * old_nnodes + os) * per + m];
				}
			delete[] *f.p;
		}
		*f.p = dst;
	}

	slot.swap(new_slot);
	col_off.swap(new_col_off);
	nkept = new_nkept;
	trash = new_nkept;
	nnodes = new_nnodes;
	compacted = true;

	const double mb_old = (double)old_nnodes / new_nnodes;
	cfdwm_mpi::out() << "lsgrid: storage compacted to " << new_nkept << " of "
		<< (long)nx * ny << " nodes (" << 100.0 * new_nkept / ((double)nx * ny)
		<< "%), " << mb_old << "x less memory -- the selected borders plus a "
		<< HALO << "-node halo (interpolation + gradient stencils, with 2 "
		"nodes of accurate margin past the retained region)."
		<< std::endl;
}

void lsGridSpline::allocate() {
	// Field arrays are sized by the STORAGE node count. Until compact()
	// runs that is the full grid, so this is the original dense layout.
	if (slot.size() != (size_t)nx * ny || nnodes == 0)
		build_identity_slot();

	// Primary fields
	ETA = new double[4 * nnodes];
	UX = new double[4 * nnodes * nl];
	UY = new double[4 * nnodes * nl];
	UZ = new double[4 * nnodes * nl];

	// secondary fields (derivatives)
	ETAdx = new double[2 * nnodes];
	ETAdy = new double[2 * nnodes];
	ETAdt = new double[4 * nnodes];
	Udt = new double[4 * nnodes * nl];
	Udx = new double[2 * nnodes * nl];
	Udy = new double[2 * nnodes * nl];
	Uds = new double[2 * nnodes * nl];
	Vdt = new double[4 * nnodes * nl];
	Vdx = new double[2 * nnodes * nl];
	Vdy = new double[2 * nnodes * nl];
	Vds = new double[2 * nnodes * nl];
	Wdt = new double[4 * nnodes * nl];
	Wdx = new double[2 * nnodes * nl];
	Wdy = new double[2 * nnodes * nl];
	Wds = new double[2 * nnodes * nl];

	// Dynamic pressure (optional: the GABC's p_nh source)
	if (store_pressure) {
		P   = new double[4 * nnodes * nl];
		Pdt = new double[4 * nnodes * nl];
		Pdx = new double[2 * nnodes * nl];
		Pdy = new double[2 * nnodes * nl];
		Pds = new double[2 * nnodes * nl];
	}

	// Misc fields
	IGNORE = new int[nnodes];
}

// A streching function for setting variable layer thickness
double lsGridSpline::slayer(int layerno) {
	//fprintf(stdout,"numlayers: %d",nl);
	double sfac = 3.0; // fixme: this should be made dimensionless and a function of specified wave
	double* dd = new double[nl];
	double ddsum = 0.;
	for (int ii = 0; ii < nl; ii++) {
		//fprintf(stdout,"%d",ii);
		dd[ii] = (1. + (nl - ii) * sfac);
		ddsum += dd[ii];
	}
	double layer_percentage = dd[layerno] / ddsum;
	delete[] dd;
	return layer_percentage;
}
template <typename T>
T clamp(const T& n, const T& lower, const T& upper) {
	return std::max(lower, std::min(n, upper));
}

// A function for constant layer thickness (equal spacing as a function of z)
double lsGridSpline::clayer(int layerno) {
	//fprintf(stdout,"numlayers: %d",nl);
	return 1./double(nl);
}

// Transforms normal z axis to streched sigma coordinates 
// s defined between -1 (seabed) and 0 (free surface)
double lsGridSpline::z2s(double z, double wave_elev, double depth) {
	return (z - (wave_elev - swl)) / (depth + (wave_elev - swl));
}

// Transforms stretched sigma coordinate to normal z
double lsGridSpline::s2z(double s, double wave_elev, double depth) {
	return (wave_elev - swl) + s * (depth + (wave_elev - swl));
}

/* Exponential vertical stretching; see layer_stretch.h for the mapping and for
   how stretch_beta is chosen when it is left automatic. */
double lsGridSpline::u2sigma(double u) {
	return layer_stretch::u2sigma(u, stretch_beta);
}

double lsGridSpline::sigma2u(double sigma) {
	return layer_stretch::sigma2u(sigma, stretch_beta);
}

double lsGridSpline::layer_thickness_ratio() const {
	return layer_stretch::thickness_ratio(stretch_beta, nl);
}

/* Automatic vertical layer distribution ('stretch_beta auto'): pick the
   stretching that minimises the vertical interpolation error for the loaded
   sea state. The mathematics lives in layer_stretch.h; this only assembles the
   component set and reports the outcome. */
void lsGridSpline::resolve_auto_stretch(Irregular& irregular) {
	if (!stretch_auto && !nl_auto) return;

	const int nc = irregular.nfreq * irregular.ndir;
	const double h = water_depth;

	// per-component horizontal velocity amplitude squared at the surface;
	// summing over direction removes any dependence on the mean heading
	std::vector<double> c, kk;
	c.reserve(std::max(nc, 0));
	kk.reserve(std::max(nc, 0));
	for (int j = 0; j < nc; j++) {
		const double amp = (irregular.omega[j] > 0.0)
			? irregular.A[j] * irregular.g * irregular.k[j] / irregular.omega[j]
			: 0.0;
		c.push_back(amp * amp);
		kk.push_back(irregular.k[j]);
	}

	/* Automatic layer count: the smallest nl whose vertical Hermite
	   reconstruction of the RMS velocity profile meets nl_tol. Under
	   stretch_auto each candidate is evaluated at its own optimal grading
	   (which then also resolves the stretching); with a prescribed beta the
	   candidates are evaluated at that grading. */
	if (nl_auto) {
		nl_auto = false;                 // resolve once
		const int NL_AUTO_MAX = 25;
		int pick = -1;
		double pick_beta = stretch_beta, pick_err = -1.0;
		for (int n = 3; n <= NL_AUTO_MAX; n++) {
			double e = -1.0, b = stretch_beta;
			if (stretch_auto) b = layer_stretch::optimal_beta(n, c, kk, h, &e);
			else e = layer_stretch::profile_error(stretch_beta, n, c, kk, h);
			if (e < 0.0 || (stretch_auto && b < 0.0)) { pick = -1; break; }
			pick = n; pick_beta = b; pick_err = e;
			if (e <= nl_tol) break;
		}
		if (pick < 0) {
			if (cfdwm_mpi::rank() == 0)
				std::cout << "lsgrid: automatic nl requested but the sea state "
					"does not define a velocity profile; keeping nl=" << nl
					<< std::endl;
			// fall through: stretch_auto (if set) resolves at the kept nl
		}
		else {
			if (pick != nl) {
				nl = pick;
				ds = 1. / std::max(1., double(nl - 1));
				/* the storage is sized by nl -> rebuild it. The IGNORE mask
				   (set_ignore ran before initialization) is nl-independent
				   and reallocate() would wipe it, so carry it across. */
				std::vector<int> keep_ignore(IGNORE, IGNORE + nnodes);
				reallocate();
				std::copy(keep_ignore.begin(), keep_ignore.end(), IGNORE);
			}
			if (cfdwm_mpi::rank() == 0)
				std::cout << "lsgrid: nl auto -> nl=" << nl
					<< " (relative L2 error of the RMS velocity profile "
					<< pick_err << ", tolerance " << nl_tol
					<< (pick_err > nl_tol ? " NOT met at the search cap" : "")
					<< ")" << std::endl;
			if (stretch_auto) {
				stretch_auto = false;   // the search already optimised it
				stretch_beta = pick_beta;
				if (cfdwm_mpi::rank() == 0)
					std::cout << "Layer stretching: auto -> beta=" << stretch_beta
						<< " (bottom/surface layer thickness ratio "
						<< layer_thickness_ratio() << " for nl=" << nl
						<< "), relative L2 error of the RMS velocity profile "
						<< pick_err << std::endl;
			}
			return;
		}
	}

	if (!stretch_auto) return;
	stretch_auto = false;                    // resolve once

	double err = 0.0;
	const double beta = layer_stretch::optimal_beta(nl, c, kk, h, &err);
	if (beta < 0.0) {
		if (cfdwm_mpi::rank() == 0)
			std::cout << "Layer stretching: auto requested but the sea state "
				"does not define a velocity profile; keeping beta="
				<< stretch_beta << std::endl;
		return;
	}
	stretch_beta = beta;

	if (cfdwm_mpi::rank() == 0) {
		std::cout << "Layer stretching: auto -> beta=" << stretch_beta
			<< " (bottom/surface layer thickness ratio "
			<< layer_thickness_ratio() << " for nl=" << nl
			<< "), relative L2 error of the RMS velocity profile "
			<< err << std::endl;
	}
}

double lsGridSpline::cart2sigmaS(double zpt, double wave_elev, double depth){
    double z = std::min(zpt, wave_elev); // impose an upper bound
    double sigma = std::max(z2s(z, wave_elev, depth), -1.); // transform to sigma coordinates
    return sigma2u(sigma);
    // defined now from 0 to 1 where 0 is sea bed and 1 is surface
}

double lsGridSpline::sigmaS2cart(double s, double wave_elev, double depth){
	// s defined now from 0 to 1 where 0 is sea bed and 1 is surface
	return s2z(u2sigma(s), wave_elev, depth);
}

void lsGridSpline::reallocate() {
	delete[] ETA; delete[] UX; delete[] UY; delete[] UZ;
	delete[] ETAdx; delete[] ETAdy; delete[] ETAdt;
	delete[] Udt; delete[] Udx; delete[] Udy; delete[] Uds;
	delete[] Vdt; delete[] Vdx; delete[] Vdy; delete[] Vds;
	delete[] Wdt; delete[] Wdx; delete[] Wdy; delete[] Wds;
	delete[] P; delete[] Pdt; delete[] Pdx; delete[] Pdy; delete[] Pds;
	P = Pdt = Pdx = Pdy = Pds = nullptr;
	delete[] IGNORE;
	allocate();
}

bool lsGridSpline::try_sigma_fft_setup(Irregular& irregular) {
	// the sigma-FFT fast path assumes a uniform depth and a spatially
	// periodic field: incompatible with variable bathymetry
	if (irregular.seabed.active()) return false;
	if (sfft.enabled) return true;
	// sigma value of each layer (x-independent by construction)
	std::vector<double> svals(nl);
	for (int m = 0; m < nl; m++) {
		svals[m] = u2sigma(ds * m);
	}
	// grid spacings; the sigma-FFT requires them to equal the lattice
	// tile node spacing (the input reader snaps the grid accordingly)
	const double dxg = (nx > 1)
		? (domain[1] - domain[0]) / (nx - 1) : 0.0;
	const double dyg = (ny > 1)
		? (domain[3] - domain[2]) / (ny - 1) : 0.0;
	return sfft.setup(irregular, dxg, dyg,
		domain[0], domain[2], svals, swl, nx, ny);
}

void lsGridSpline::sigma_fft_fill_level(Irregular& irregular, double t,
                                        int tt) {
	std::vector<double> eta_v, u_v, v_v, w_v, p_v;
	sfft.fill(irregular, t, eta_v, u_v, v_v, w_v,
		store_pressure ? &p_v : nullptr);
	const int nfx = sfft.nfftx, nfy = sfft.nffty;
	const int NP = nfx * nfy;
	for (int i = 0; i < nx; i++) {
		const int ii = i % nfx;       // node nx-1 wraps onto node 0
		for (int j = 0; j < ny; j++) {
			if (!fill_node(i, j)) continue;   // not stored / inside the ignore box
			const int jj = sfft.twodim ? (j % nfy) : 0;
			const int p = ii * nfy + jj;
			ETA[tt * nnodes + nid(i, j)] = eta_v[p] + swl;
			for (int m = 0; m < nl; m++) {
				const long idx = ((long)tt * nnodes + nid(i, j)) * nl + m;
				UX[idx] = u_v[(size_t)m * NP + p];
				UY[idx] = v_v[(size_t)m * NP + p];
				UZ[idx] = w_v[(size_t)m * NP + p];
				if (store_pressure) P[idx] = p_v[(size_t)m * NP + p];
			}
		}
	}
}

void lsGridSpline::update_gradient_dt(double* data, double* graddt){
    // function for computation of temporal gradients using central difference
    if (tstep > 0){
        int tid = (tstep + 2) % 4;
        int tid1 = (tstep + 3) % 4;
        int tid0 = (tstep + 1) % 4;
		#pragma omp for
		for (int i = 0; i < nnodes * nl; i++) {
            graddt[tid * nnodes * nl + i] = (data[tid1 * nnodes * nl + i] - data[tid0 * nnodes * nl + i]) / (2 * dt);
		}
	}
    else {
		#pragma omp for
        for (int i = 0; i < nnodes * nl; i++) {
            graddt[1 * nnodes * nl + i] = (data[2 * nnodes * nl + i] - data[1 * nnodes * nl + i]) / dt;
            graddt[2 * nnodes * nl + i] = (data[3 * nnodes * nl + i] - data[1 * nnodes * nl + i]) / (2 * dt);
		}
	}
}

void lsGridSpline::update_gradient_eta_dt(double* data, double* graddt){
    // function for computation of temporal gradients using central difference
    if (tstep > 0){
        int tid = (tstep + 2) % 4;
        int tid1 = (tstep + 3) % 4;
        int tid0 = (tstep + 1) % 4;

		#pragma omp for
        for (int i = 0; i < nnodes; i++) {
            graddt[tid * nnodes + i] = (data[tid1 * nnodes + i] - data[tid0 * nnodes + i]) / (2 * dt);
		}
	}
    else {
		#pragma omp for
        for (int i = 0; i < nnodes; i++) {
            graddt[1 * nnodes + i] = (data[2 * nnodes + i] - data[1 * nnodes + i]) / dt;
            graddt[2 * nnodes + i] = (data[3 * nnodes + i] - data[1 * nnodes + i]) / (2 * dt);
		}
	}
}

void lsGridSpline::update_gradient_eta_dxdy(double* eta, double* gradx,double* grady){
    // tid: specifies which time step to compute gradients for (value from 0 to 1)
    // compute spatial gradient of the sea surface using central difference

    int tid0 = (tstep + 1) % 4;
    int tid1 = (tstep + 2) % 4;

	if (nx == 1){
		#pragma omp for
		for (int j = 0; j < ny; j++) {
            gradx[0 * nnodes + nid(0, j)] = 0.;
            gradx[1 * nnodes + nid(0, j)] = 0.;
		}
	}
	else{
		#pragma omp for
		for (int i = 1; i < (nx-1); i++) {
			for (int j = 0; j < ny; j++) {
				gradx[0 * nnodes + nid(i, j)] = (eta[tid0 * nnodes + nid(i + 1, j)] -
													eta[tid0 * nnodes + nid(i - 1, j)]) / (2 * dx);
				gradx[1 * nnodes + nid(i, j)] = (eta[tid1 * nnodes + nid(i + 1, j)] -
													eta[tid1 * nnodes + nid(i - 1, j)]) / (2 * dx);
			}
		}
		
		// boundaries are special cases (x=0 and x=xn)
		#pragma omp for
		for (int j = 0; j < ny; j++) {
			gradx[0 * nnodes + nid(0, j)] = (eta[tid0 * nnodes + nid(1, j)] - eta[tid0 * nnodes + nid(0, j)]) / dx;
			gradx[0 * nnodes + nid(nx-1, j)] = (eta[tid0 * nnodes + nid(nx-1, j)] -
													eta[tid0 * nnodes + nid(nx-2, j)]) / dx;
			gradx[1 * nnodes + nid(0, j)] = (eta[tid1 * nnodes + nid(1, j)] - eta[tid1 * nnodes + nid(0, j)]) / dx;
			gradx[1 * nnodes + nid(nx - 1, j)] = (eta[tid1 * nnodes + nid(nx - 1, j)] -
													eta[tid1 * nnodes + nid(nx - 2, j)]) / dx;
		}
	}

	if (ny == 1){
		#pragma omp for
		for (int i = 0; i < nx; i++) {
            grady[0 * nnodes + nid(i, 0)] = 0.;
            grady[1 * nnodes + nid(i, 0)] = 0.;
		}
	}
	else{
		// compute spatial gradient using central difference
		#pragma omp for
		for (int i = 0; i < nx; i++) {
			for (int j = 1; j < (ny-1); j++) {
				grady[0 * nnodes + nid(i, j)] = (eta[tid0 * nnodes + nid(i, j + 1)] -
													eta[tid0 * nnodes + nid(i, j - 1)]) / (2 * dy);
				grady[1 * nnodes + nid(i, j)] = (eta[tid1 * nnodes + nid(i, j + 1)] -
													eta[tid1 * nnodes + nid(i, j - 1)]) / (2 * dy);
			}
		}
		#pragma omp for
		// y=0 and y=yn
		for (int i = 0; i < nx; i++) {
			grady[0 * nnodes + nid(i, 0)] = (eta[tid0 * nnodes + nid(i, 1)] - eta[tid0 * nnodes + nid(i, 0)]) / dy;
			grady[0 * nnodes + nid(i, ny - 1)] = (eta[tid0 * nnodes + nid(i, ny - 1)] -
														eta[tid0 * nnodes + nid(i, ny - 2)]) / dy;
			grady[1 * nnodes + nid(i, 0)] = (eta[tid1 * nnodes + nid(i, 1)] - eta[tid1 * nnodes + nid(i, 0)]) / dy;
			grady[1 * nnodes + nid(i, ny - 1)] = (eta[tid1 * nnodes + nid(i, ny - 1)] -
													eta[tid1 * nnodes + nid(i, ny - 2)]) / dy;

		}
	}
}

void lsGridSpline::update_gradient_dxdydz(double* data, double* gradx, double* grady, double* gradz){
    // compute spatial gradient using central difference

    int tid0 = (tstep + 1) % 4;
    int tid1 = (tstep + 2) % 4;
	if (nx == 1){
		#pragma omp for
		for (int j = 0; j < ny; j++) {
			for (int m = 0; m < nl; m++) {
				gradx[(0 * nnodes + nid(0, j)) * nl + m] = 0.;
				gradx[(1 * nnodes + nid(0, j)) * nl + m] = 0.;
			}
		}
	}
	else{
		#pragma omp for collapse(2)
		for (int i = 1; i < (nx-1); i++) {
			for (int j = 0; j < ny; j++) {
				for (int m = 0; m < nl; m++) {
					gradx[(0 * nnodes + nid(i, j)) * nl + m] = ((data[(tid0 * nnodes + nid(i + 1, j)) * nl + m] -
						data[(tid0 * nnodes + nid(i - 1, j)) * nl + m]) / (2 * dx));
					gradx[(1 * nnodes + nid(i, j)) * nl + m] = ((data[(tid1 * nnodes + nid(i + 1, j)) * nl + m] -
						data[(tid1 * nnodes + nid(i - 1, j)) * nl + m]) / (2 * dx));
				}
			}
		}
		// x=0 and x=xn
		#pragma omp for collapse(2)
		for (int j = 0; j < ny; j++) {
			for (int m = 0; m < nl; m++) {
				// nid(0, j), NOT the dense j*nl: after compaction the
				// storage index of column 0 is no longer the grid index,
				// and the dense form writes out of (or across) the
				// compacted arrays
				gradx[(0 * nnodes + nid(0, j)) * nl + m] = (
					data[(tid0 * nnodes + nid(1, j)) * nl + m] -
					data[(tid0 * nnodes + nid(0, j)) * nl + m]) / dx;
				gradx[(0 * nnodes + nid(nx-1, j)) * nl + m] = (
					data[(tid0 * nnodes + nid(nx-1, j)) * nl + m] -
					data[(tid0 * nnodes + nid(nx-2, j)) * nl + m]) / dx;
				gradx[(1 * nnodes + nid(0, j)) * nl + m] = (
					data[(tid1 * nnodes + nid(1, j)) * nl + m] -
					data[(tid1 * nnodes + nid(0, j)) * nl + m]) / dx;
				gradx[(1 * nnodes + nid(nx - 1, j)) * nl + m] = (
					data[(tid1 * nnodes + nid(nx - 1, j)) * nl + m] -
					data[(tid1 * nnodes + nid(nx - 2, j)) * nl + m]) / dx;
			}
		}
	}
	if (ny == 1){
		#pragma omp for
		for (int i = 0; i < nx; i++) {
			for (int m = 0; m < nl; m++) {
				// nid(i, 0), NOT the dense i*nl: the dense form runs off
				// the end of the compacted arrays (heap corruption on any
				// 1D border-mode grid)
				grady[(0 * nnodes + nid(i, 0)) * nl + m] = 0.;
				grady[(1 * nnodes + nid(i, 0)) * nl + m] = 0.;
			}
		}
	}
	else{
		// compute spatial gradient using central difference
		#pragma omp for collapse(2)
		for (int i = 0; i < nx; i++) {
			for (int j = 1; j < (ny-1); j++) {
				for (int m = 0; m < nl; m++) {
					grady[(0 * nnodes + nid(i, j)) * nl + m] = (
							(data[(tid0 * nnodes + nid(i, j + 1)) * nl + m] -
							data[(tid0 * nnodes + nid(i, j - 1)) * nl + m]) / (2 * dy));
					grady[(1 * nnodes + nid(i, j)) * nl + m] = (
							(data[(tid1 * nnodes + nid(i, j + 1)) * nl + m] -
							data[(tid1 * nnodes + nid(i, j - 1)) * nl + m]) / (2 * dy));
				}
			}
		}

		// y=0 and y=yn
		#pragma omp for collapse(2)
		for (int i = 0; i < nx; i++) {
			for (int m = 0; m < nl; m++) {
				grady[(0 * nnodes + nid(i, 0)) * nl + m] = (
						(data[(tid0 * nnodes + nid(i, 1)) * nl + m] -
						data[(tid0 * nnodes + nid(i, 0)) * nl + m]) / dy);
				grady[(0 * nnodes + nid(i, ny - 1)) * nl + m] = (
						(data[(tid0 * nnodes + nid(i, ny - 1)) * nl + m] -
						data[(tid0 * nnodes + nid(i, ny - 2)) * nl + m]) / dy);
				grady[(1 * nnodes + nid(i, 0)) * nl + m] = (
						(data[(tid1 * nnodes + nid(i, 1)) * nl + m] -
						data[(tid1 * nnodes + nid(i, 0)) * nl + m]) / dy);
				grady[(1 * nnodes + nid(i, ny - 1)) * nl + m] = (
						(data[(tid1 * nnodes + nid(i, ny - 1)) * nl + m] -
						data[(tid1 * nnodes + nid(i, ny - 2)) * nl + m]) / dy);
			}
		}
	}
    if (nl == 1){
		#pragma omp for
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				gradz[0 * nnodes + nid(i, j)] = 0.;
				gradz[1 * nnodes + nid(i, j)] = 0.;
			}
		}
	}
	else{
		// compute grad z spatial gradient using central difference
		#pragma omp for collapse(2)
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				for (int m = 1; m < (nl-1); m++) {
					gradz[(0 * nnodes + nid(i, j)) * nl + m] = (
							(data[(tid0 * nnodes + nid(i, j)) * nl + (m + 1)] -
							data[(tid0 * nnodes + nid(i, j)) * nl + (m - 1)]) / (2 * ds));
					gradz[(1 * nnodes + nid(i, j)) * nl + m] = (
							(data[(tid1 * nnodes + nid(i, j)) * nl + (m + 1)] -
							data[(tid1 * nnodes + nid(i, j)) * nl + (m - 1)]) / (2 * ds));
				}
			}
		}
		// z=0 and z=zn
		#pragma omp for collapse(2)
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				gradz[(0 * nnodes + nid(i, j)) * nl] = (
						(data[(tid0 * nnodes + nid(i, j)) * nl + 1] -
						data[(tid0 * nnodes + nid(i, j)) * nl]) / ds);
				gradz[(0 * nnodes + nid(i, j)) * nl + (nl - 1)] = (
						(data[(tid0 * nnodes + nid(i, j)) * nl + (nl - 1)] -
						data[(tid0 * nnodes + nid(i, j)) * nl + (nl - 2)]) / ds);
				gradz[(1 * nnodes + nid(i, j)) * nl] = (
						(data[(tid1 * nnodes + nid(i, j)) * nl + 1] -
						data[(tid1 * nnodes + nid(i, j)) * nl]) / ds);
				gradz[(1 * nnodes + nid(i, j)) * nl + (nl - 1)] = (
						(data[(tid1 * nnodes + nid(i, j)) * nl + (nl - 1)] -
						data[(tid1 * nnodes + nid(i, j)) * nl + (nl - 2)]) / ds);
			}
		}
	}
}

void lsGridSpline::square_vals(double* C, double* data, int nxp, int nyp, int tid){
    const int ip = clamp(nxp + 1, 0, nx - 1), jp = clamp(nyp + 1, 0, ny - 1);
    C[0] = data[tid * nnodes + nid(nxp, nyp)]; // C00
    C[1] = data[tid * nnodes + nid(nxp, jp)];  // C01
    C[2] = data[tid * nnodes + nid(ip,  nyp)]; // C10
    C[3] = data[tid * nnodes + nid(ip,  jp)];  // C11
}

void lsGridSpline::cube_vals(double* C, double* data, int nxp, int nyp, int nlp, int tid){
    const int ip = clamp(nxp + 1, 0, nx - 1), jp = clamp(nyp + 1, 0, ny - 1);
    const int mp = clamp(nlp + 1, 0, nl - 1);
    const long a = (long)tid * nnodes + nid(nxp, nyp);
    const long b = (long)tid * nnodes + nid(nxp, jp);
    const long c = (long)tid * nnodes + nid(ip,  nyp);
    const long d = (long)tid * nnodes + nid(ip,  jp);
    C[0] = data[a * nl + nlp]; //C000
    C[1] = data[b * nl + nlp]; //C010
    C[2] = data[c * nl + nlp]; //C100
    C[3] = data[d * nl + nlp]; //C110
    C[4] = data[a * nl + mp];  //C001
    C[5] = data[b * nl + mp];  //C011
    C[6] = data[c * nl + mp];  //C101
    C[7] = data[d * nl + mp];  //C111
}

double lsGridSpline::spline_interp_velo(double* U, double* Udt, double* Udx, double* Udy, double* Uds, int nxp, int nyp, int nsp0, int nsp1, double xd, double yd, double sd0, double sd1, double td, int tid1, int tid2){
    // Ux
    // --------------------------------------------------------------------
    // u0, u1, u0dt, u1dt, u0dx, u1dx, u0dy, u1dy, u0ds, u1ds
	double C[8], D[8], Ct[8], Dt[8], Cx[8], Dx[8], Cy[8], Dy[8], Cs[8], Ds[8];
    cube_vals(C, U, nxp, nyp, nsp0, tid1);
    cube_vals(D, U, nxp, nyp, nsp1, tid2);
    cube_vals(Ct, Udt, nxp, nyp, nsp0, tid1);
    cube_vals(Dt, Udt, nxp, nyp, nsp1, tid2);
    cube_vals(Cx, Udx, nxp, nyp, nsp0, 0);
    cube_vals(Dx, Udx, nxp, nyp, nsp1, 1);
    cube_vals(Cy, Udy, nxp, nyp, nsp0, 0);
    cube_vals(Dy, Udy, nxp, nyp, nsp1, 1);
    cube_vals(Cs, Uds, nxp, nyp, nsp0, 0);
    cube_vals(Ds, Uds, nxp, nyp, nsp1, 1);

    // gcompute spline coefficients
    double a00 = dx * Cx[0] - (C[2] - C[0]);
    double b00 = -dx * Cx[2] + (C[2] - C[0]);
    double a10 = dx * Cx[1] - (C[3] - C[1]);
    double b10 = -dx * Cx[3] + (C[3] - C[1]);
    double a01 = dx * Cx[4] - (C[6] - C[4]);
    double b01 = -dx * Cx[6] + (C[6] - C[4]);
	double a11 = dx * Cx[5] - (C[7] - C[5]);
    double b11 = -dx * Cx[7] + (C[7] - C[5]);
    double c00 = dx * Dx[0] - (D[2] - D[0]);
    double d00 = -dx * Dx[2] + (D[2] - D[0]);
    double c10 = dx * Dx[1] - (D[3] - D[1]);
    double d10 = -dx * Dx[3] + (D[3] - D[1]);
    double c01 = dx * Dx[4] - (D[6] - D[4]);
    double d01 = -dx * Dx[6] + (D[6] - D[4]);
	double c11 = dx * Dx[5] - (D[7] - D[5]);
    double d11 = -dx * Dx[7] + (D[7] - D[5]);

    // compute new values
    double C00 = C[0] * (1. - xd) + C[2] * xd + xd * (1 - xd) * (a00 * (1 - xd) + b00 * xd);
    double C10 = C[1] * (1. - xd) + C[3] * xd + xd * (1 - xd) * (a10 * (1 - xd) + b10 * xd);
    double C01 = C[4] * (1. - xd) + C[6] * xd + xd * (1 - xd) * (a01 * (1 - xd) + b01 * xd);
    double C11 = C[5] * (1. - xd) + C[7] * xd + xd * (1 - xd) * (a11 * (1 - xd) + b11 * xd);
    double D00 = D[0] * (1. - xd) + D[2] * xd + xd * (1 - xd) * (c00 * (1 - xd) + d00 * xd);
    double D10 = D[1] * (1. - xd) + D[3] * xd + xd * (1 - xd) * (c10 * (1 - xd) + d10 * xd);
    double D01 = D[4] * (1. - xd) + D[6] * xd + xd * (1 - xd) * (c01 * (1 - xd) + d01 * xd);
    double D11 = D[5] * (1. - xd) + D[7] * xd + xd * (1 - xd) * (c11 * (1 - xd) + d11 * xd);

    // Reduce y, s and t
    double C00y = Cy[0] * (1. - xd) + Cy[2] * xd;
    double C10y = Cy[1] * (1. - xd) + Cy[3] * xd;
    double C01y = Cy[4] * (1. - xd) + Cy[6] * xd;
    double C11y = Cy[5] * (1. - xd) + Cy[7] * xd;
    double D00y = Dy[0] * (1. - xd) + Dy[2] * xd;
    double D10y = Dy[1] * (1. - xd) + Dy[3] * xd;
    double D01y = Dy[4] * (1. - xd) + Dy[6] * xd;
    double D11y = Dy[5] * (1. - xd) + Dy[7] * xd;

    double C00t = Ct[0] * (1. - xd) + Ct[2] * xd;
    double C10t = Ct[1] * (1. - xd) + Ct[3] * xd;
    double C01t = Ct[4] * (1. - xd) + Ct[6] * xd;
    double C11t = Ct[5] * (1. - xd) + Ct[7] * xd;
    double D00t = Dt[0] * (1. - xd) + Dt[2] * xd;
    double D10t = Dt[1] * (1. - xd) + Dt[3] * xd;
    double D01t = Dt[4] * (1. - xd) + Dt[6] * xd;
    double D11t = Dt[5] * (1. - xd) + Dt[7] * xd;

    double C00s = Cs[0] * (1. - xd) + Cs[2] * xd;
    double C10s = Cs[1] * (1. - xd) + Cs[3] * xd;
    double C01s = Cs[4] * (1. - xd) + Cs[6] * xd;
    double C11s = Cs[5] * (1. - xd) + Cs[7] * xd;
    double D00s = Ds[0] * (1. - xd) + Ds[2] * xd;
    double D10s = Ds[1] * (1. - xd) + Ds[3] * xd;
    double D01s = Ds[4] * (1. - xd) + Ds[6] * xd;
    double D11s = Ds[5] * (1. - xd) + Ds[7] * xd;

    // compute spline coefficients
    double a0 = dy * C00y - (C10 - C00);
    double b0 = -dy * C10y + (C10 - C00);
    double a1 = dy * C01y - (C11 - C01);
    double b1 = -dy * C11y + (C11 - C01);
    double c0 = dy * D00y - (D10 - D00);
    double d0 = -dy * D10y + (D10 - D00);
    double c1 = dy * D01y - (D11 - D01);
    double d1 = -dy * D11y + (D11 - D01);

    // compute new values
    double C0 = C00 * (1. - yd) + C10 * yd + yd * (1 - yd) * (a0 * (1 - yd) + b0 * yd);
    double C1 = C01 * (1. - yd) + C11 * yd + yd * (1 - yd) * (a1 * (1 - yd) + b1 * yd);
    double D0 = D00 * (1. - yd) + D10 * yd + yd * (1 - yd) * (c0 * (1 - yd) + d0 * yd);
    double D1 = D01 * (1. - yd) + D11 * yd + yd * (1 - yd) * (c1 * (1 - yd) + d1 * yd);

    // reduce s and t direction
    double C0s = C00s * (1. - yd) + C10s * yd;
    double C1s = C01s * (1. - yd) + C11s * yd;
    double D0s = D00s * (1. - yd) + D10s * yd;
    double D1s = D01s * (1. - yd) + D11s * yd;
    double C0t = C00t * (1. - yd) + C10t * yd;
    double C1t = C01t * (1. - yd) + C11t * yd;
    double D0t = D00t * (1. - yd) + D10t * yd;
    double D1t = D01t * (1. - yd) + D11t * yd;

    // compute spline coefficients
    a0 = ds * C0s - (C1 - C0);
    b0 = -ds * C1s + (C1 - C0);
    a1 = ds * D0s - (D1 - D0);
    b1 = -ds * D1s + (D1 - D0);

    // compute wave elevation for two time steps
    double CC = C0 * (1. - sd0) + C1 * sd0 + sd0 * (1 - sd0) * (a0 * (1 - sd0) + b0 * sd0);
    double DD = D0 * (1. - sd1) + D1 * sd1 + sd1 * (1 - sd1) * (a1 * (1 - sd1) + b1 * sd1);

    // compute spline coefficients
    double CCt = C0t * (1. - sd0) + C1t * sd0;
    double DDt = D0t * (1. - sd1) + D1t * sd1;
    double a = dt * CCt - (DD - CC);
    double b = -dt * DDt + (DD - CC);

    // compute final interpolated wave elevation for single point
    return CC * (1. - td) + DD * td + td * (1 - td) * (a * (1 - td) + b * td);
}


std::vector<double> lsGridSpline::get_kinematics_at_point(double tpt, double xpt, double ypt, double zpt, double h, bool with_pressure){
    // for t, for x, for y, for z stacking

    // keep the caller's coordinates: the wrap below rewrites xpt/ypt, and an
    // error reported at the wrapped position is misleading -- it is not where
    // the caller asked
    const double xq0 = xpt, yq0 = ypt;

    // periodic wrap: the equal-dk field repeats with the domain lengths,
    // so any point maps back into the domain
    if (periodic_x) {
        const double Lp = domain[1] - domain[0];
        double xr = std::fmod(xpt - domain[0], Lp);
        if (xr < 0.0) xr += Lp;
        xpt = domain[0] + xr;
    }
    if (periodic_y) {
        const double Lp = domain[3] - domain[2];
        double yr = std::fmod(ypt - domain[2], Lp);
        if (yr < 0.0) yr += Lp;
        ypt = domain[2] + yr;
    }

    // list of required data
    // eta0, eta1, eta1dt, eta2dt, eta1dx, eta2dx, eta1dy, eta2dy
    // u0, u1, u0dt, u1dt, u0dx, u1dx, u0dy, u1dy, u0ds, u1ds
    // v0, v1, v0dt, v1dt, v0dx, v1dx, v0dy, v1dy, v0ds, v1ds
    // w0, w1, w0dt, w1dt, w0dx, w1dx, w0dy, w1dy, w0ds, w1ds

    int tid1 = (tstep + 1) % 4;
    int tid2 = (tstep + 2) % 4;

    // Step 1 - Interpolate to find wave elevation at point (xpt,ypt)
	float nxp_temp, nyp_temp;
	double xd = std::modf(clamp((xpt - domain[0]) / dx, 0., nx - 1.), &nxp_temp);
	double yd = std::modf(clamp((ypt - domain[2]) / dy, 0., ny - 1.), &nyp_temp);
	double td = std::min(1., std::max(0., (tpt - t0) / dt));
	
	int nxp = int(nxp_temp);
	int nyp = int(nyp_temp);

	// After compaction the nodes under the ignore box hold no data at all.
	// Reading them would silently return the scratch node, so refuse: a
	// query here means the ignore box does not match what the solver asks
	// for. (Before compaction every node is stored and this never trips.)
	if (compacted) {
		const int ip = clamp(nxp + 1, 0, nx - 1), jp = clamp(nyp + 1, 0, ny - 1);
		// grad_ok, not merely stored: the Hermite blend reads the PRECOMPUTED
		// gradients at each corner, and a corner on the outermost kept ring
		// has neighbours that were dropped, so its gradient is meaningless.
		// Checking storage alone would let a query two nodes inside the box
		// through and quietly interpolate with a garbage slope.
		if (!grad_ok(nxp, nyp) || !grad_ok(nxp, jp) ||
		    !grad_ok(ip, nyp) || !grad_ok(ip, jp)) {
			// name what IS stored, so the mismatch is actionable.
			// NB rank-divergent failure: std::cerr directly (per the
			// mpi_utils.h rule) -- cfdwm_mpi::err() is muted on
			// non-root ranks and this fires on whichever rank owns
			// the querying subdomain, silently killing the job.
			std::cerr << "[rank " << cfdwm_mpi::rank()
				<< "] lsgrid ERROR: kinematics requested at ("
				<< xq0 << ", " << yq0 << ")"
				<< ((xq0 != xpt || yq0 != ypt)
				    ? std::string(", which the periodic tile wraps to (")
				      + std::to_string(xpt) + ", " + std::to_string(ypt) + ")"
				    : std::string())
				<< ", which is not stored. Kinematics "
				"are kept only along the selected borders ("
				<< (border[0] ? "left " : "") << (border[1] ? "right " : "")
				<< (border[2] ? "bottom " : "") << (border[3] ? "top " : "")
				<< "), " << BORDER_NODES << " nodes deep plus the stencil "
				"halo. Either enable the border this point belongs to, or "
				"query only near the borders that are enabled." << std::endl;
			cfdwm_mpi::abort(1);
		}
	}

	double C[4], D[4], Ct[4], Dt[4], Cx[4], Dx[4], Cy[4], Dy[4], Cs[4], Ds[4];
    square_vals(C, ETA, nxp, nyp, tid1);
    square_vals(D, ETA, nxp, nyp, tid2);
    square_vals(Ct, ETAdt, nxp, nyp, tid1);
    square_vals(Dt, ETAdt, nxp, nyp, tid2);
    square_vals(Cx, ETAdx, nxp, nyp, 0);
    square_vals(Dx, ETAdx, nxp, nyp, 1);
    square_vals(Cy, ETAdy, nxp, nyp, 0);
    square_vals(Dy, ETAdy, nxp, nyp, 1);

	// compute spline coefficients
	double a00 = dx * Cx[0] - (C[2] - C[0]);
	double a01 = dx * Cx[1] - (C[3] - C[1]);
	double b00 = -dx * Cx[2] + (C[2] - C[0]);
	double b01 = -dx * Cx[3] + (C[3] - C[1]);
	double a10 = dx * Dx[0] - (D[2] - D[0]);
	double a11 = dx * Dx[1] - (D[3] - D[1]);
	double b10 = -dx * Dx[2] + (D[2] - D[0]);
	double b11 = -dx * Dx[3] + (D[3] - D[1]);

	// compute new values
	double C0 = C[0] * (1. - xd) + C[2] * xd + xd * (1 - xd) * (a00 * (1 - xd) + b00 * xd);
	double C1 = C[1] * (1. - xd) + C[3] * xd + xd * (1 - xd) * (a01 * (1 - xd) + b01 * xd);
	double D0 = D[0] * (1. - xd) + D[2] * xd + xd * (1 - xd) * (a10 * (1 - xd) + b10 * xd);
	double D1 = D[1] * (1. - xd) + D[3] * xd + xd * (1 - xd) * (a11 * (1 - xd) + b11 * xd);

    // reduce y and t direction
    double C0y = Cy[0] * (1. - xd) + Cy[2] * xd;
    double C1y = Cy[1] * (1. - xd) + Cy[3] * xd;
    double D0y = Dy[0] * (1. - xd) + Dy[2] * xd;
    double D1y = Dy[1] * (1. - xd) + Dy[3] * xd;
    double C0t = Ct[0] * (1. - xd) + Ct[2] * xd;
    double C1t = Ct[1] * (1. - xd) + Ct[3] * xd;
    double D0t = Dt[0] * (1. - xd) + Dt[2] * xd;
    double D1t = Dt[1] * (1. - xd) + Dt[3] * xd;

    // compute spline coefficients
    double a0 = dy * C0y - (C1 - C0);
    double b0 = -dy * C1y + (C1 - C0);
    double a1 = dy * D0y - (D1 - D0);
    double b1 = -dy * D1y + (D1 - D0);

    // compute wave elevation for two time steps
    double wave_elev0 = C0 * (1. - yd) + C1 * yd + yd * (1 - yd) * (a0 * (1 - yd) + b0 * yd);
    double wave_elev1 = D0 * (1. - yd) + D1 * yd + yd * (1 - yd) * (a1 * (1 - yd) + b1 * yd);

    // compute spline coefficients
    double CCt = C0t * (1. - yd) + C1t * yd;
    double DDt = D0t * (1. - yd) + D1t * yd;
    double a = dt * CCt - (wave_elev1 - wave_elev0);
    double b = -dt * DDt + (wave_elev1 - wave_elev0);

    // compute final interpolated wave elevation for single point
    double wave_elev = wave_elev0 * (1. - td) + wave_elev1 * td + td * (1 - td) * (a * (1 - td) + b * td);

    // Step 2 - Now compute vertical sigma coords and find kinematics from SigmaStreched grid.

    double spt0 = cart2sigmaS(zpt, wave_elev0, local_depth(xpt, ypt));
    double nsp0a = spt0 / ds;
    double spt1 = cart2sigmaS(zpt, wave_elev1, local_depth(xpt, ypt));
    double nsp1a = spt1 / ds;

	float sd0_temp, sd1_temp;
    double sd0 = std::modf(nsp0a, &sd0_temp);
    double sd1 = std::modf(nsp1a, &sd1_temp);
    int nsp0 = int(sd0_temp);
    int nsp1 = int(sd1_temp);

    double ux = spline_interp_velo(UX, Udt, Udx, Udy, Uds, nxp, nyp, nsp0, nsp1, xd, yd, sd0, sd1, td, tid1, tid2);
    double vy = spline_interp_velo(UY, Vdt, Vdx, Vdy, Vds, nxp, nyp, nsp0, nsp1, xd, yd, sd0, sd1, td, tid1, tid2);
    double wz = spline_interp_velo(UZ, Wdt, Wdx, Wdy, Wds, nxp, nyp, nsp0, nsp1, xd, yd, sd0, sd1, td, tid1, tid2);

	std::vector<double> res;
	res.push_back(wave_elev);
	res.push_back(ux);
	res.push_back(vy);
	res.push_back(wz);
	// only meaningful when the grid carries a pressure field; the caller
	// (wave_nhPres / wave_DynPres) checks store_pressure before asking
	if (with_pressure && store_pressure)
		res.push_back(spline_interp_velo(P, Pdt, Pdx, Pdy, Pds, nxp, nyp, nsp0, nsp1, xd, yd, sd0, sd1, td, tid1, tid2));

	return res;
}

//#define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))


bool lsGridSpline::CheckTime(double tpt) {
	/* Checks to see if the time tpt is within the interval t0 to t1. If so, returns true*/
	if (tpt > t0 + dt) {
		if (cfdwm_mpi::rank() == 0)
			std::cout << "t0: " << t0 << ", t1: " << (t0 + dt) << ", tpt: " << tpt << std::endl;
		return true;
	}
	return false;
}

// function to find if given point 
// lies inside a given rectangle or not. 
bool lsGridSpline::CheckBounds()
{
	if (periodic_x && periodic_y) return true;   // fully periodic
	if (periodic_x) {   // x always in bounds by periodic wrap
		if (bymin >= domain[2] && bymax <= domain[3]) return true;
		std::cout << "Requested point outside specified grid domain (y). "
			"adjust the bounds of the grid and try again." << std::endl;
		cfdwm_mpi::abort(-1);
	}
	if (bxmin >= domain[0] && bxmax <= domain[1] && bymin >= domain[2] && bymax <= domain[3])
		return true;
	else {
		std::cout << "Requested point outside specified grid domain. adjust the bounds of the grid and try again." << std::endl;
		cfdwm_mpi::abort(-1);
		return false;
	}
}

void lsGridSpline::update_bounds(double xpt, double ypt) {
	bxmin = std::min(xpt, bxmin);
	bxmax = std::max(xpt, bxmax);
	bymin = std::min(ypt, bymin);
	bymax = std::max(ypt, bymax);

}


//----------------------------------------------------------------------------------------------------------------------------------------


// Second order theory functions


//----------------------------------------------------------------------------------------------------------------------------------------

// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 3
// MPI: work-balanced column partition (by active-cell count). Computed once after
// the IGNORE mask is set; reused by initialize_kinematics and update.
void lsGridSpline::compute_mpi_partition() {
	int s = cfdwm_mpi::size();
	if (s <= 1) { mpi_bounds.assign(2, 0); mpi_bounds[1] = nx; mpi_i0 = 0; mpi_i1 = nx; return; }
	std::vector<long> colw(nx);
	for (int i = 0; i < nx; i++) {
		long a = 0;
		for (int j = 0; j < ny; j++) a += (IGNORE[i * ny + j] ? 0 : 1);
		colw[i] = a;
	}
	cfdwm_mpi::partition_weighted(colw, nx, mpi_bounds);
	int r = cfdwm_mpi::rank();
	mpi_i0 = mpi_bounds[r]; mpi_i1 = mpi_bounds[r + 1];
}

void lsGridSpline::initialize_kinematics(Irregular& irregular) {
	fill_all = true;

	bxmin = domain[0];
	bxmax = domain[1];
	bymin = domain[2];
	bymax = domain[3];

	dx = (domain[1] - domain[0]) / std::max(1., double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));
	//dz = (domain_end[2] - domain_start[2]) / double(NZ - 1);
	ds = 1. / std::max(1., double(nl - 1));

	// must precede anything that consumes the sigma mapping
	resolve_auto_stretch(irregular);

	double tt0[4] = {t0-dt, t0, t0+dt, t0+2*dt};

	double dd = omp_get_wtime();

	// sigma-FFT fast path (equal-dk + smit): whole time levels via FFT
	if (try_sigma_fft_setup(irregular)) {
		compute_mpi_partition();   // queried later; fill is rank-redundant
		for (int tt = 0; tt < 4; tt++) {
			sigma_fft_fill_level(irregular, tt0[tt], tt);
		}
		update_gradient_eta_dxdy(ETA, ETAdx, ETAdy);
		update_gradient_eta_dt(ETA, ETAdt);
		update_gradient_dxdydz(UX, Udx, Udy, Uds);
		update_gradient_dt(UX, Udt);
		update_gradient_dxdydz(UY, Vdx, Vdy, Vds);
		update_gradient_dt(UY, Vdt);
		update_gradient_dxdydz(UZ, Wdx, Wdy, Wds);
		update_gradient_dt(UZ, Wdt);
		if (store_pressure) {
			update_gradient_dxdydz(P, Pdx, Pdy, Pds);
			update_gradient_dt(P, Pdt);
		}
		if (dump_vtk) write_vtk(false);
		if (cfdwm_mpi::rank() == 0)
			std::cout << "Generation of domain kinematics data using the "
				"sigma-FFT layer fill completed. Initialization time: "
				<< (omp_get_wtime() - dd) << " seconds." << std::endl;
		return;
	}

	//omp_set_num_threads(1)
	//omp_set_num_threads(omp_get_max_threads());

	compute_mpi_partition();   // MPI: work-balanced column partition
	                           // (before the parallel region: it resizes
	                           // shared vectors and must run exactly once)
#pragma omp parallel // start parallel initialization
	{
#pragma omp master
		if (cfdwm_mpi::rank() == 0)
			std::cout << "Number of available threads: " << omp_get_num_threads() << std::endl;

		// Main grid
		int i0 = mpi_i0, i1 = mpi_i1;
#pragma omp for collapse(2)
		for (int tt = 0; tt < 4; tt++){ // We initialize 4 steps in time
			for (int i = i0; i < i1; i++) {
				for (int j = 0; j < ny; j++) {
					double xpt = domain[0] + dx * i;
					double ypt = domain[2] + dy * j;

					if (irregular.seabed.active()) {
						// variable bathymetry: slope-corrected first order
						// + depth-local S&D second order; layers follow
						// the local depth
						ETA[tt * nnodes + nid(i, j)] = irregular.eta1_slope(tt0[tt], xpt, ypt) + swl;
						if (irregular.order >= 2)
							ETA[tt * nnodes + nid(i, j)] += irregular.eta2_slope(tt0[tt], xpt, ypt);
						double eta_temp = ETA[tt * nnodes + nid(i, j)];
						for (int m = 0; m < nl; m++) {
							double zpt0 = sigmaS2cart(ds*m, eta_temp, local_depth(xpt, ypt));
							UVWP S = irregular.uvwp1_slope(tt0[tt], xpt, ypt, zpt0 - swl);
							if (irregular.order >= 2) {
								UVWP S2 = irregular.uvwp2_slope(tt0[tt], xpt, ypt, zpt0 - swl);
								S.u += S2.u; S.v += S2.v; S.w += S2.w; S.p += S2.p;
							}
							UX[(tt * nnodes + nid(i, j)) * nl + m] = S.u;
							UY[(tt * nnodes + nid(i, j)) * nl + m] = S.v;
							UZ[(tt * nnodes + nid(i, j)) * nl + m] = S.w;
							if (store_pressure)
								P[(tt * nnodes + nid(i, j)) * nl + m] = S.p;
						}
						continue;
					}

					ETA[tt * nnodes + nid(i, j)] = irregular.eta1(tt0[tt], xpt, ypt) + irregular.eta2(tt0[tt], xpt, ypt) + swl;

					double eta_temp = ETA[tt * nnodes + nid(i, j)];

					double PHI0_dxdz = irregular.phi1_dxdz(tt0[tt], xpt, ypt);
					double PHI0_dydz = irregular.phi1_dydz(tt0[tt], xpt, ypt);
					double PHI0_dzdz = irregular.phi1_dzdz(tt0[tt], xpt, ypt);

					for (int m = 0; m < nl; m++) {
						double zpt0 = sigmaS2cart(ds*m, eta_temp, local_depth(xpt, ypt));
						double zpt0max = std::max(0., zpt0 - swl);
						double zpt0min = std::min(0., zpt0); // swl already included in irregular class functions

						//UX0[i * ny * nl + j * nl + m] = irregular.u1(t0, xpt, ypt, zpt0min) + irregular.u2(t0, xpt, ypt, zpt0min) + PHI0_dxdz * zpt0max;
						//UY0[i * ny * nl + j * nl + m] = irregular.v1(t0, xpt, ypt, zpt0min) + irregular.v2(t0, xpt, ypt, zpt0min) + PHI0_dydz * zpt0max;
						//UZ0[i * ny * nl + j * nl + m] = irregular.w1(t0, xpt, ypt, zpt0min) + irregular.w2(t0, xpt, ypt, zpt0min) + PHI0_dzdz * zpt0max;

						UVWP U2 = irregular.uvw2(tt0[tt], xpt, ypt, zpt0min);
						UX[(tt * nnodes + nid(i, j)) * nl + m] = irregular.u1(tt0[tt], xpt, ypt, zpt0min) + U2.u + PHI0_dxdz * zpt0max;
						UY[(tt * nnodes + nid(i, j)) * nl + m] = irregular.v1(tt0[tt], xpt, ypt, zpt0min) + U2.v + PHI0_dydz * zpt0max;
						UZ[(tt * nnodes + nid(i, j)) * nl + m] = irregular.w1(tt0[tt], xpt, ypt, zpt0min) + U2.w + PHI0_dzdz * zpt0max;
						if (store_pressure)
							P[(tt * nnodes + nid(i, j)) * nl + m] = irregular.dp(tt0[tt], xpt, ypt, zpt0);
					}
				}
			}
		}
	}

	// MPI: assemble the full grid from each rank's column block, then compute
	// gradients (redundant on every rank; cheap finite differences).
	for (int tt = 0; tt < 4; tt++) {
		cfdwm_mpi::allgather_columns_var(ETA + (long)tt * nnodes, mpi_bounds, col_off, (long)1);
		cfdwm_mpi::allgather_columns_var(UX  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
		cfdwm_mpi::allgather_columns_var(UY  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
		cfdwm_mpi::allgather_columns_var(UZ  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
		if (store_pressure)
			cfdwm_mpi::allgather_columns_var(P + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
	}
	// Update derivatives
	update_gradient_eta_dxdy(ETA, ETAdx, ETAdy);
	update_gradient_eta_dt(ETA, ETAdt);
	update_gradient_dxdydz(UX, Udx, Udy, Uds);
	update_gradient_dt(UX, Udt);
	update_gradient_dxdydz(UY, Vdx, Vdy, Vds);
	update_gradient_dt(UY, Vdt);
	update_gradient_dxdydz(UZ, Wdx, Wdy, Wds);
	update_gradient_dt(UZ, Wdt);
	if (store_pressure) {
		update_gradient_dxdydz(P, Pdx, Pdy, Pds);
		update_gradient_dt(P, Pdt);
	}

	// End parallel initialization
	if (dump_vtk) {
		write_vtk(false);
	}
	dd = omp_get_wtime() - dd;
	if (cfdwm_mpi::rank() == 0) {
		std::cout << "Generation of domain kinematics data using irregular second order wave theory completed. ";
		std::cout << "Initialization time: " << dd << " seconds." << std::endl;
		std::cout << "Interpolation can commence..." << std::endl;
	}
}

void lsGridSpline::initialize_kinematics_with_ignore(Irregular& irregular) {
	fill_all = false;

	bxmin = domain[0];
	bxmax = domain[1];
	bymin = domain[2];
	bymax = domain[3];

	dx = (domain[1] - domain[0]) / std::max(1., double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));
	//dz = (domain_end[2] - domain_start[2]) / double(NZ - 1);
	ds = 1. / std::max(1., double(nl - 1));

	// must precede anything that consumes the sigma mapping
	resolve_auto_stretch(irregular);

	double tt0[4] = {t0, t0, t0+dt, t0+2*dt};

	double dd = omp_get_wtime();

	// sigma-FFT fast path: filling everything is cheaper than skipping
	if (try_sigma_fft_setup(irregular)) {
		// Say so explicitly. An FFT produces the whole lattice in one
		// transform, so restricting storage to the borders cannot skip any
		// of that work -- it only saves memory here (it saves work on the
		// direct-summation path). Reading the log without this, the border
		// selection looks silently disregarded.
		long nign = 0;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				if (IGNORE[i * ny + j]) nign++;
		const long ntot = (long)nx * ny;   // grid nodes, not storage slots
		if (nign)
			cfdwm_mpi::out()
			    << "sigma-fft: " << (ntot - nign) << " of " << ntot
			    << " lsgrid nodes are inside the selected borders. NOTE: the "
			       "sigma-FFT computes the whole lattice in a single "
			       "transform, so all " << ntot << " nodes are evaluated "
			       "regardless -- restricting storage to the borders saves "
			       "memory, not FFT work (it saves work on the "
			       "direct-summation path)." << std::endl;
		compute_mpi_partition();
		for (int tt = 0; tt < 4; tt++) {
			sigma_fft_fill_level(irregular, tt0[tt], tt);
		}
		update_gradient_eta_dxdy(ETA, ETAdx, ETAdy);
		update_gradient_eta_dt(ETA, ETAdt);
		update_gradient_dxdydz(UX, Udx, Udy, Uds);
		update_gradient_dt(UX, Udt);
		update_gradient_dxdydz(UY, Vdx, Vdy, Vds);
		update_gradient_dt(UY, Vdt);
		update_gradient_dxdydz(UZ, Wdx, Wdy, Wds);
		update_gradient_dt(UZ, Wdt);
		if (store_pressure) {
			update_gradient_dxdydz(P, Pdx, Pdy, Pds);
			update_gradient_dt(P, Pdt);
		}
		if (dump_vtk) write_vtk(false);
		if (cfdwm_mpi::rank() == 0)
			std::cout << "Generation of domain kinematics data using the "
				"sigma-FFT layer fill completed. Initialization time: "
				<< (omp_get_wtime() - dd) << " seconds." << std::endl;
		return;
	}

	//omp_set_num_threads(1)
	//omp_set_num_threads(omp_get_max_threads());

	compute_mpi_partition();   // MPI: work-balanced column partition
	                           // (before the parallel region: it resizes
	                           // shared vectors and must run exactly once)
#pragma omp parallel // start parallel initialization
	{
#pragma omp master
		if (cfdwm_mpi::rank() == 0)
			std::cout << "Number of available threads: " << omp_get_num_threads() << std::endl;

		// Main grid
		int i0 = mpi_i0, i1 = mpi_i1;
#pragma omp for collapse(2)
		for (int tt = 0; tt < 4; tt++){ // We initialize 4 steps in time
			for (int i = i0; i < i1; i++) {
				for (int j = 0; j < ny; j++) {
					if (fill_node(i, j)) {
						double xpt = domain[0] + dx * i;
						double ypt = domain[2] + dy * j;

						if (irregular.seabed.active()) {
							// variable bathymetry: slope-corrected first order
							// + depth-local S&D second order; layers
							// follow the local depth
							ETA[tt * nnodes + nid(i, j)] = irregular.eta1_slope(tt0[tt], xpt, ypt) + swl;
							if (irregular.order >= 2)
								ETA[tt * nnodes + nid(i, j)] += irregular.eta2_slope(tt0[tt], xpt, ypt);
							double eta_temp = ETA[tt * nnodes + nid(i, j)];
							for (int m = 0; m < nl; m++) {
								double zpt0 = sigmaS2cart(ds*m, eta_temp, local_depth(xpt, ypt));
								UVWP S = irregular.uvwp1_slope(tt0[tt], xpt, ypt, zpt0 - swl);
								if (irregular.order >= 2) {
									UVWP S2 = irregular.uvwp2_slope(tt0[tt], xpt, ypt, zpt0 - swl);
									S.u += S2.u; S.v += S2.v; S.w += S2.w; S.p += S2.p;
								}
								UX[(tt * nnodes + nid(i, j)) * nl + m] = S.u;
								UY[(tt * nnodes + nid(i, j)) * nl + m] = S.v;
								UZ[(tt * nnodes + nid(i, j)) * nl + m] = S.w;
								if (store_pressure)
									P[(tt * nnodes + nid(i, j)) * nl + m] = S.p;
							}
							continue;
						}

						ETA[tt * nnodes + nid(i, j)] = irregular.eta1(tt0[tt], xpt, ypt) + irregular.eta2(tt0[tt], xpt, ypt) + swl;

						double eta_temp = ETA[tt * nnodes + nid(i, j)];

						double PHI0_dxdz = irregular.phi1_dxdz(tt0[tt], xpt, ypt);
						double PHI0_dydz = irregular.phi1_dydz(tt0[tt], xpt, ypt);
						double PHI0_dzdz = irregular.phi1_dzdz(tt0[tt], xpt, ypt);

						for (int m = 0; m < nl; m++) {
							double zpt0 = sigmaS2cart(ds*m, eta_temp, local_depth(xpt, ypt));
							double zpt0max = std::max(0., zpt0 - swl);
							double zpt0min = std::min(0., zpt0); // swl already included in irregular class functions

							//UX0[i * ny * nl + j * nl + m] = irregular.u1(t0, xpt, ypt, zpt0min) + irregular.u2(t0, xpt, ypt, zpt0min) + PHI0_dxdz * zpt0max;
							//UY0[i * ny * nl + j * nl + m] = irregular.v1(t0, xpt, ypt, zpt0min) + irregular.v2(t0, xpt, ypt, zpt0min) + PHI0_dydz * zpt0max;
							//UZ0[i * ny * nl + j * nl + m] = irregular.w1(t0, xpt, ypt, zpt0min) + irregular.w2(t0, xpt, ypt, zpt0min) + PHI0_dzdz * zpt0max;

							UVWP U2 = irregular.uvw2(tt0[tt], xpt, ypt, zpt0min);
							UX[(tt * nnodes + nid(i, j)) * nl + m] = irregular.u1(tt0[tt], xpt, ypt, zpt0min) + U2.u + PHI0_dxdz * zpt0max;
							UY[(tt * nnodes + nid(i, j)) * nl + m] = irregular.v1(tt0[tt], xpt, ypt, zpt0min) + U2.v + PHI0_dydz * zpt0max;
							UZ[(tt * nnodes + nid(i, j)) * nl + m] = irregular.w1(tt0[tt], xpt, ypt, zpt0min) + U2.w + PHI0_dzdz * zpt0max;
							if (store_pressure)
								P[(tt * nnodes + nid(i, j)) * nl + m] = irregular.dp(tt0[tt], xpt, ypt, zpt0);
						}
					}
				}
			}
		}
	}

	// MPI: assemble the full grid from each rank's column block, then compute
	// gradients (redundant on every rank; cheap finite differences).
	for (int tt = 0; tt < 4; tt++) {
		cfdwm_mpi::allgather_columns_var(ETA + (long)tt * nnodes, mpi_bounds, col_off, (long)1);
		cfdwm_mpi::allgather_columns_var(UX  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
		cfdwm_mpi::allgather_columns_var(UY  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
		cfdwm_mpi::allgather_columns_var(UZ  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
		if (store_pressure)
			cfdwm_mpi::allgather_columns_var(P + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
	}
	// Update derivatives
	update_gradient_eta_dxdy(ETA, ETAdx, ETAdy);
	update_gradient_eta_dt(ETA, ETAdt);
	update_gradient_dxdydz(UX, Udx, Udy, Uds);
	update_gradient_dt(UX, Udt);
	update_gradient_dxdydz(UY, Vdx, Vdy, Vds);
	update_gradient_dt(UY, Vdt);
	update_gradient_dxdydz(UZ, Wdx, Wdy, Wds);
	update_gradient_dt(UZ, Wdt);
	if (store_pressure) {
		update_gradient_dxdydz(P, Pdx, Pdy, Pds);
		update_gradient_dt(P, Pdt);
	}

	// End parallel initialization
	if (dump_vtk) {
		write_vtk(false);
	}
	dd = omp_get_wtime() - dd;
	if (cfdwm_mpi::rank() == 0) {
		std::cout << "Generation of domain kinematics data using irregular second order wave theory completed. ";
		std::cout << "Initialization time: " << dd << " seconds." << std::endl;
		std::cout << "Interpolation can commence..." << std::endl;
	}
}

// When called, updates the arrays storing surface elevation and kinematics data for timestep t0 = t1, t1 = t1+dt
void lsGridSpline::update(Irregular& irregular, double t_target)
{
	fill_all = false;
	// ignore_at_init == false: the whole grid was filled at t = 0 for the
	// benefit of a solver initialising its full domain. That has now served
	// its purpose -- from here on only the region outside the ignore box is
	// updated, so the rest of the storage can go back.
	if (!compacted) compact();

	double dd, t1, t2;
	int tt;

	// Start by checking bounds
	/*
	if (!disable_checkbounds){
		CheckBounds();
	}*/
	
	if (cfdwm_mpi::rank() == 0)
		std::cout << "t_target: " << t_target << "sec" << std::endl;
	tstep++;
	// new time ste
	t0 += dt; // base step
	// Updating surface elevations
		
	t1 = t0 + dt; 
	t2 = t0 + 2*dt;
	tt = (tstep + 3) % 4; // index to write new time step (two steps ahead of basestep)
	dd = omp_get_wtime();

	// sigma-FFT fast path: fill the new time slice with FFTs
	if (sfft.enabled) {
		sigma_fft_fill_level(irregular, t2, tt);
#pragma omp parallel
		{
			update_gradient_eta_dxdy(ETA, ETAdx, ETAdy);
			update_gradient_eta_dt(ETA, ETAdt);
			update_gradient_dxdydz(UX, Udx, Udy, Uds);
			update_gradient_dt(UX, Udt);
			update_gradient_dxdydz(UY, Vdx, Vdy, Vds);
			update_gradient_dt(UY, Vdt);
			update_gradient_dxdydz(UZ, Wdx, Wdy, Wds);
			update_gradient_dt(UZ, Wdt);
			if (store_pressure) {
				update_gradient_dxdydz(P, Pdx, Pdy, Pds);
				update_gradient_dt(P, Pdt);
			}
		}
		if (dump_vtk) write_vtk(false);
		dd = omp_get_wtime() - dd;
		if (cfdwm_mpi::rank() == 0) {
			std::cout << "update time: " << dd << " sec" << std::endl;
			std::cout << "lsGridSpline matrices updated. t = " << t0 << " to "
				<< (t0 + dt) << std::endl;
		}
		return;
	}

	//omp_set_num_threads(1);
	//omp_set_num_threads(omp_get_max_threads());
	int i0 = mpi_i0, i1 = mpi_i1;               // MPI: reuse the partition from init
	#pragma omp parallel // start parallel initialization
	{
		// Main grid
#pragma omp for collapse(2)
		for (int j = 0; j < ny; j++) {
			for (int i = i0; i < i1; i++) {
				//std::cout << "j: " << j << std::endl;
				double xpt = domain[0] + dx * i;
				//std::cout << "processornum: " << omp_get_thread_num() << std::endl;
				double ypt = domain[2] + dy * j;
				if (fill_node(i, j)) {
					if (irregular.seabed.active()) {
						// variable bathymetry: slope-corrected first order
						// + depth-local S&D second order; layers
						// follow the local depth
						ETA[tt * nnodes + nid(i, j)] = irregular.eta1_slope(t2, xpt, ypt) + swl;
						if (irregular.order >= 2)
							ETA[tt * nnodes + nid(i, j)] += irregular.eta2_slope(t2, xpt, ypt);
						double eta_temp = ETA[tt * nnodes + nid(i, j)];
						for (int m = 0; m < nl; m++) {
							double zpt0 = sigmaS2cart(ds*m, eta_temp, local_depth(xpt, ypt));
							UVWP S = irregular.uvwp1_slope(t2, xpt, ypt, zpt0 - swl);
							if (irregular.order >= 2) {
								UVWP S2 = irregular.uvwp2_slope(t2, xpt, ypt, zpt0 - swl);
								S.u += S2.u; S.v += S2.v; S.w += S2.w; S.p += S2.p;
							}
							UX[(tt * nnodes + nid(i, j)) * nl + m] = S.u;
							UY[(tt * nnodes + nid(i, j)) * nl + m] = S.v;
							UZ[(tt * nnodes + nid(i, j)) * nl + m] = S.w;
							if (store_pressure)
								P[(tt * nnodes + nid(i, j)) * nl + m] = S.p;
						}
						continue;
					}

					ETA[tt * nnodes + nid(i, j)] = irregular.eta1(t2, xpt, ypt) + irregular.eta2(t2, xpt, ypt) + swl;
					double eta_temp = ETA[tt * nnodes + nid(i, j)];

					double PHI1_dxdz = irregular.phi1_dxdz(t2, xpt, ypt);
					double PHI1_dydz = irregular.phi1_dydz(t2, xpt, ypt);
					double PHI1_dzdz = irregular.phi1_dzdz(t2, xpt, ypt);

					for (int m = 0; m < nl; m++) {
						double zpt0 = sigmaS2cart(ds*m, eta_temp, local_depth(xpt, ypt));
						double zpt1max = std::max(0., zpt0 - swl);
						double zpt1min = std::min(0., zpt0); // swl already included in irregular class functions

						UVWP U2 = irregular.uvw2(t2, xpt, ypt, zpt1min);
						UX[(tt * nnodes + nid(i, j)) * nl + m] = irregular.u1(t2, xpt, ypt, zpt1min) + U2.u + PHI1_dxdz * zpt1max;
						UY[(tt * nnodes + nid(i, j)) * nl + m] = irregular.v1(t2, xpt, ypt, zpt1min) + U2.v + PHI1_dydz * zpt1max;
						UZ[(tt * nnodes + nid(i, j)) * nl + m] = irregular.w1(t2, xpt, ypt, zpt1min) + U2.w + PHI1_dzdz * zpt1max;
						if (store_pressure)
							P[(tt * nnodes + nid(i, j)) * nl + m] = irregular.dp(t2, xpt, ypt, zpt0);
					}
				}
			}
		}
		
	}   // close the compute parallel region

	// MPI: assemble the new time slice from each rank's column block, then
	// recompute gradients (redundant on every rank; cheap finite differences).
	cfdwm_mpi::allgather_columns_var(ETA + (long)tt * nnodes, mpi_bounds, col_off, (long)1);
	cfdwm_mpi::allgather_columns_var(UX  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
	cfdwm_mpi::allgather_columns_var(UY  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
	cfdwm_mpi::allgather_columns_var(UZ  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
	if (store_pressure)
		cfdwm_mpi::allgather_columns_var(P + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);

	#pragma omp parallel
	{
		update_gradient_eta_dxdy(ETA, ETAdx, ETAdy);
		update_gradient_eta_dt(ETA, ETAdt);
		update_gradient_dxdydz(UX, Udx, Udy, Uds);
		update_gradient_dt(UX, Udt);
		update_gradient_dxdydz(UY, Vdx, Vdy, Vds);
		update_gradient_dt(UY, Vdt);
		update_gradient_dxdydz(UZ, Wdx, Wdy, Wds);
		update_gradient_dt(UZ, Wdt);   // was missing: Wdt went stale after updates
		if (store_pressure) {
			update_gradient_dxdydz(P, Pdx, Pdy, Pds);
			update_gradient_dt(P, Pdt);
		}
	}

	if (dump_vtk) {
		write_vtk(false);
	}
	dd = omp_get_wtime() - dd;
	if (cfdwm_mpi::rank() == 0) {
		std::cout << "update time: " << dd << " sec" << std::endl;
		std::cout << "lsGridSpline matrices updated. t = " << t0 << " to " << (t0 + dt) << std::endl;
	}

}

//----------------------------------------------------------------------------------------------------------------------------------------


// SWD Routines


//----------------------------------------------------------------------------------------------------------------------------------------


#if defined(SWD_enable)

// SWD as kinematics source.
//
// The SWD library is stateful: UpdateTime() mutates per-object interpolation
// state and reads spectral timesteps from the file, and the bundled build
// supports only ONE open swd object per process (fixed Fortran file unit --
// a concurrent second open fails and can even corrupt the live object's file
// stream). All swd sampling below is therefore strictly SERIAL within the
// process: one UpdateTime per time level, then a plain double loop over the
// columns. MPI ranks are separate processes with their own swd objects, so
// the column partition still parallelizes the fill across ranks.
//
// Sampling: eta and the grid layers are stored in the CFD frame (swl
// included); the swd object is queried in the SWD frame (z=0 at the calm
// surface) and the sample point is clamped at the instantaneous surface,
// since SWD kinematics diverge exponentially above z=eta.

void lsGridSpline::initialize_kinematics(SpectralWaveData* swd) {

	bxmin = domain[0];
	bxmax = domain[1];
	bymin = domain[2];
	bymax = domain[3];

	dx = (domain[1] - domain[0]) / std::max(1., double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));
	ds = 1. / std::max(1., double(nl - 1));

	// Level 0 sits one step before the start time, but swd time cannot be
	// negative: clamp at 0 (duplicated first level -- same startup
	// approximation as the with_ignore variant).
	double tt0[4] = {std::max(t0 - dt, 0.), t0, t0 + dt, t0 + 2 * dt};

	double dd = omp_get_wtime();

	compute_mpi_partition();                     // MPI: work-balanced column partition
	int i0 = mpi_i0, i1 = mpi_i1;

	for (int tt = 0; tt < 4; tt++) { // We initialize 4 steps in time
		try {
			swd->UpdateTime(tt0[tt]);
		}
		catch (SwdException& e) {  // Could be t > tmax from file.
			std::cerr << typeid(e).name() << std::endl << e.what() << std::endl;
			cfdwm_mpi::abort(EXIT_FAILURE);
		}
		for (int i = i0; i < i1; i++) {
			for (int j = 0; j < ny; j++) {
				double xpt = domain[0] + dx * i;
				double ypt = domain[2] + dy * j;

				double eta_swd = swd->Elev(xpt, ypt);   // SWD frame (z=0 at calm surface)
				ETA[tt * nnodes + nid(i, j)] = eta_swd + swl;

				double eta_temp = ETA[tt * nnodes + nid(i, j)];

				for (int m = 0; m < nl; m++) {
					double zpt0 = sigmaS2cart(ds * m, eta_temp, local_depth(xpt, ypt)); // CFD frame
					double zpt_swd = std::min(zpt0 - swl, eta_swd);           // SWD frame, clamped at surface
					vector_swd U = swd->GradPhi(xpt, ypt, zpt_swd);
					UX[(tt * nnodes + nid(i, j)) * nl + m] = U.x;
					UY[(tt * nnodes + nid(i, j)) * nl + m] = U.y;
					UZ[(tt * nnodes + nid(i, j)) * nl + m] = U.z;
				}
			}
		}
	}

	// MPI: assemble the full grid from each rank's column block, then compute
	// gradients (redundant on every rank; cheap finite differences).
	for (int tt = 0; tt < 4; tt++) {
		cfdwm_mpi::allgather_columns_var(ETA + (long)tt * nnodes, mpi_bounds, col_off, (long)1);
		cfdwm_mpi::allgather_columns_var(UX  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
		cfdwm_mpi::allgather_columns_var(UY  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
		cfdwm_mpi::allgather_columns_var(UZ  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
	}
	// Update derivatives (no swd access in here)
	update_gradient_eta_dxdy(ETA, ETAdx, ETAdy);
	update_gradient_eta_dt(ETA, ETAdt);
	update_gradient_dxdydz(UX, Udx, Udy, Uds);
	update_gradient_dt(UX, Udt);
	update_gradient_dxdydz(UY, Vdx, Vdy, Vds);
	update_gradient_dt(UY, Vdt);
	update_gradient_dxdydz(UZ, Wdx, Wdy, Wds);
	update_gradient_dt(UZ, Wdt);

	if (dump_vtk) {
		write_vtk(false);
	}
	dd = omp_get_wtime() - dd;
	if (cfdwm_mpi::rank() == 0) {
		std::cout << "Generation of domain kinematics data using SWD completed. ";
		std::cout << "Initialization time: " << dd << " seconds." << std::endl;
		std::cout << "Interpolation can commence..." << std::endl;
	}
}

void lsGridSpline::initialize_kinematics_with_ignore(SpectralWaveData* swd) {

	bxmin = domain[0];
	bxmax = domain[1];
	bymin = domain[2];
	bymax = domain[3];

	dx = (domain[1] - domain[0]) / std::max(1., double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));
	ds = 1. / std::max(1., double(nl - 1));

	double tt0[4] = {t0, t0, t0 + dt, t0 + 2 * dt};

	double dd = omp_get_wtime();

	compute_mpi_partition();                     // MPI: work-balanced column partition
	int i0 = mpi_i0, i1 = mpi_i1;

	for (int tt = 0; tt < 4; tt++) { // We initialize 4 steps in time
		try {
			swd->UpdateTime(tt0[tt]);
		}
		catch (SwdException& e) {  // Could be t > tmax from file.
			std::cerr << typeid(e).name() << std::endl << e.what() << std::endl;
			cfdwm_mpi::abort(EXIT_FAILURE);
		}
		for (int i = i0; i < i1; i++) {
			for (int j = 0; j < ny; j++) {
				if (fill_node(i, j)) {
					double xpt = domain[0] + dx * i;
					double ypt = domain[2] + dy * j;

					double eta_swd = swd->Elev(xpt, ypt);   // SWD frame (z=0 at calm surface)
					ETA[tt * nnodes + nid(i, j)] = eta_swd + swl;

					double eta_temp = ETA[tt * nnodes + nid(i, j)];

					for (int m = 0; m < nl; m++) {
						double zpt0 = sigmaS2cart(ds * m, eta_temp, local_depth(xpt, ypt)); // CFD frame
						double zpt_swd = std::min(zpt0 - swl, eta_swd);           // SWD frame, clamped at surface
						vector_swd U = swd->GradPhi(xpt, ypt, zpt_swd);
						UX[(tt * nnodes + nid(i, j)) * nl + m] = U.x;
						UY[(tt * nnodes + nid(i, j)) * nl + m] = U.y;
						UZ[(tt * nnodes + nid(i, j)) * nl + m] = U.z;
					}
				}
			}
		}
	}

	// MPI: assemble the full grid from each rank's column block, then compute
	// gradients (redundant on every rank; cheap finite differences).
	for (int tt = 0; tt < 4; tt++) {
		cfdwm_mpi::allgather_columns_var(ETA + (long)tt * nnodes, mpi_bounds, col_off, (long)1);
		cfdwm_mpi::allgather_columns_var(UX  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
		cfdwm_mpi::allgather_columns_var(UY  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
		cfdwm_mpi::allgather_columns_var(UZ  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
	}
	// Update derivatives (no swd access in here)
	update_gradient_eta_dxdy(ETA, ETAdx, ETAdy);
	update_gradient_eta_dt(ETA, ETAdt);
	update_gradient_dxdydz(UX, Udx, Udy, Uds);
	update_gradient_dt(UX, Udt);
	update_gradient_dxdydz(UY, Vdx, Vdy, Vds);
	update_gradient_dt(UY, Vdt);
	update_gradient_dxdydz(UZ, Wdx, Wdy, Wds);
	update_gradient_dt(UZ, Wdt);

	if (dump_vtk) {
		write_vtk(false);
	}
	dd = omp_get_wtime() - dd;
	if (cfdwm_mpi::rank() == 0) {
		std::cout << "Generation of domain kinematics data using SWD completed. ";
		std::cout << "Initialization time: " << dd << " seconds." << std::endl;
		std::cout << "Interpolation can commence..." << std::endl;
	}
}

// When called, updates the arrays storing surface elevation and kinematics
// data for timestep t0 = t1, t1 = t1+dt (swd variant; serial fill, see above)
void lsGridSpline::update(SpectralWaveData* swd, double t_target)
{
	double dd, t2;
	int tt;

	if (cfdwm_mpi::rank() == 0)
		std::cout << "t_target: " << t_target << "sec" << std::endl;
	tstep++;
	// new time step
	t0 += dt; // base step
	t2 = t0 + 2 * dt;
	tt = (tstep + 3) % 4; // index to write new time step (two steps ahead of basestep)
	dd = omp_get_wtime();
	int i0 = mpi_i0, i1 = mpi_i1;               // MPI: reuse the partition from init

	try {
		swd->UpdateTime(t2);
	}
	catch (SwdException& e) {  // Could be t > tmax from file.
		std::cerr << typeid(e).name() << std::endl << e.what() << std::endl;
		cfdwm_mpi::abort(EXIT_FAILURE);
	}
	for (int i = i0; i < i1; i++) {
		for (int j = 0; j < ny; j++) {
			if (fill_node(i, j)) {
				double xpt = domain[0] + dx * i;
				double ypt = domain[2] + dy * j;

				double eta_swd = swd->Elev(xpt, ypt);   // SWD frame (z=0 at calm surface)
				ETA[tt * nnodes + nid(i, j)] = eta_swd + swl;

				double eta_temp = ETA[tt * nnodes + nid(i, j)];

				for (int m = 0; m < nl; m++) {
					double zpt0 = sigmaS2cart(ds * m, eta_temp, local_depth(xpt, ypt)); // CFD frame
					double zpt_swd = std::min(zpt0 - swl, eta_swd);           // SWD frame, clamped at surface
					vector_swd U = swd->GradPhi(xpt, ypt, zpt_swd);
					UX[(tt * nnodes + nid(i, j)) * nl + m] = U.x;
					UY[(tt * nnodes + nid(i, j)) * nl + m] = U.y;
					UZ[(tt * nnodes + nid(i, j)) * nl + m] = U.z;
				}
			}
		}
	}

	// MPI: assemble the new time slice from each rank's column block, then
	// recompute gradients (redundant on every rank; cheap finite differences).
	cfdwm_mpi::allgather_columns_var(ETA + (long)tt * nnodes, mpi_bounds, col_off, (long)1);
	cfdwm_mpi::allgather_columns_var(UX  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
	cfdwm_mpi::allgather_columns_var(UY  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);
	cfdwm_mpi::allgather_columns_var(UZ  + (long)tt * nnodes * nl, mpi_bounds, col_off, (long)nl);

	#pragma omp parallel
	{
		update_gradient_eta_dxdy(ETA, ETAdx, ETAdy);
		update_gradient_eta_dt(ETA, ETAdt);
		update_gradient_dxdydz(UX, Udx, Udy, Uds);
		update_gradient_dt(UX, Udt);
		update_gradient_dxdydz(UY, Vdx, Vdy, Vds);
		update_gradient_dt(UY, Vdt);
		update_gradient_dxdydz(UZ, Wdx, Wdy, Wds);
		update_gradient_dt(UZ, Wdt);
	}

	if (dump_vtk) {
		write_vtk(false);
	}
	dd = omp_get_wtime() - dd;
	if (cfdwm_mpi::rank() == 0) {
		std::cout << "update time: " << dd << " sec" << std::endl;
		std::cout << "lsGridSpline matrices updated. t = " << t0 << " to " << (t0 + dt) << std::endl;
	}
}

#endif // SWD_enable

//----------------------------------------------------------------------------------------------------------------------------------------


// VTK output


//----------------------------------------------------------------------------------------------------------------------------------------


// when called, writes stored kinematics to file
void lsGridSpline::write_vtk(bool endtime) {
	// MPI: every rank holds the full assembled grid; one writer suffices
	// (all ranks writing the same file concurrently would corrupt it)
	if (cfdwm_mpi::rank() != 0) return;
	char buffer[256]; sprintf(buffer, "%05d", tstep);

	if (dirExists(vtk_directory_path.c_str()) == 0) {
		std::cout << "WARNING: Specified directory for storage of VTK files does not exist. Directory will be created at the following path:  " << vtk_directory_path << std::endl;
		createDirectory(vtk_directory_path);
	}


	std::string str(buffer);
	std::string fpath = (vtk_directory_path + vtk_prefix + buffer + ".vts");
	std::cout << fpath << std::endl;
	FILE* fp = fopen(fpath.c_str(), "w");
	if (endtime)
		export_vts(fp, true);
	else
		export_vts(fp, false);
	fclose(fp);

	std::cout << "wrote kinematics to: " << fpath << std::endl;
}

/* exports sGrid at t = t0 (or t0 + dt) to a .vts structured-grid file for
   visualization in vtk/paraview. The lsgrid is logically structured
   (nx x ny x nl, sigma-stretched z), so a StructuredGrid carries the
   geometry without connectivity tables. VTK point ordering: x varies
   fastest, then y, then the layer index (z). */
void lsGridSpline::export_vts(FILE* fp, bool last)
{
	int tt = (tstep + 1) % 4; // index for base step
	double tval = t0;
	if (last) {
		tt = (tstep + 2) % 4; // index for basestep + 1
		tval = t0 + dt;
	}

	// write header
	fputs("<?xml version=\"1.0\"?>\n"
		"<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
	fprintf(fp, "\t <StructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n", nx - 1, ny - 1, nl - 1);
	fprintf(fp, "\t\t <FieldData> \n");
	fprintf(fp, "\t\t\t <DataArray type = \"Float64\" Name = \"%s\" NumberOfTuples = \"1\" format = \"ascii\" RangeMin = \"%.3f\" RangeMax = \"%.3f\"> \n", vtk_timelabel.c_str(), tval, tval);
	fprintf(fp, "\t\t\t %.3f \n", tval);
	fprintf(fp, "\t\t\t </DataArray > \n");
	fprintf(fp, "\t\t </FieldData> \n");

	fprintf(fp, "\t\t <Piece Extent=\"0 %d 0 %d 0 %d\">\n", nx - 1, ny - 1, nl - 1);

	fputs("\t\t\t <PointData Scalars=\"scalars\">\n", fp);
	fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n");
	for (int m = 0; m < nl; m++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				fprintf(fp, "%g %g %g\n", UX[(tt * nnodes + nid(i, j)) * nl + m], UY[(tt * nnodes + nid(i, j)) * nl + m], UZ[(tt * nnodes + nid(i, j)) * nl + m]);
			}
		}
	}
	fputs("\t\t\t\t </DataArray>\n", fp);
	if (store_pressure) {
		fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"1\" Name=\"PressureDyn\" format=\"ascii\">\n");
		for (int m = 0; m < nl; m++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
					fprintf(fp, "%g\n", P[(tt * nnodes + nid(i, j)) * nl + m]);
				}
			}
		}
		fputs("\t\t\t\t </DataArray>\n", fp);
	}
	fputs("\t\t\t </PointData>\n", fp);

	fputs("\t\t\t <Points>\n", fp);
	fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
	for (int m = 0; m < nl; m++) {
		double upt = ds * m;
		for (int j = 0; j < ny; j++) {
			double ypt = domain[2] + dy * j;
			for (int i = 0; i < nx; i++) {
				double xpt = domain[0] + dx * i;
				double eta_temp = ETA[tt * nnodes + nid(i, j)];
				double zpt = s2z(u2sigma(upt), eta_temp, local_depth(xpt, ypt));
				fprintf(fp, "%12.4f %12.4f %12.4f\n", xpt, ypt, zpt);
			}
		}
	}
	fputs("\t\t\t\t </DataArray>\n", fp);
	fputs("\t\t\t </Points>\n", fp);

	fputs("\t\t </Piece>\n", fp);
	fputs("\t </StructuredGrid>\n", fp);
	fputs("</VTKFile>\n", fp);
	fflush(fp);
}

/* Set area of domain to ignore when update kinematics data. this is useful when prescribing kinematics at the boundaries*/
void lsGridSpline::set_ignore()
{
	dx = (domain[1] - domain[0]) / std::max(1.,double(nx - 1));
	dy = (domain[3] - domain[2]) / std::max(1., double(ny - 1));

	// QA: the numbers the grid actually ended up with, whether they came
	// from the input file or were derived. dx/dy/dt are what matter for
	// resolution and are not otherwise visible -- only nx/ny/dt are stated.
	cfdwm_mpi::out() << "lsgrid: " << nx << " x " << ny << " x " << nl
		<< " nodes over [" << domain[0] << ", " << domain[1] << "] x ["
		<< domain[2] << ", " << domain[3] << "] -- dx=" << dx << " m, dy="
		<< dy << " m, dt=" << dt << " s." << std::endl;

	// Keep a strip along each selected border, drop everything else.
	// Everything downstream -- compact(), fill_node(), the query guard --
	// consumes this mask, so borders are expressed here and nowhere else.
	// With no border selected nothing is dropped and the grid stays dense,
	// which is what a client initialising its whole domain wants.
	// No border selected (border mode off, or on with all four zero) means
	// store the whole grid -- what a client initialising its full domain
	// wants. Dropping everything instead leaves an empty grid.
	if (!border_mode || !(border[0] || border[1] || border[2] || border[3])) {
		for (long k = 0; k < (long)nx * ny; k++) IGNORE[k] = 0;
		return;
	}
	// Size the strip in METRES when the client has declared how far inward
	// it queries: sizing in nodes alone silently stops being enough once the
	// client's cells are coarse relative to dx.
	int nb = BORDER_NODES;
	if (border_margin > 0.) {
		const double dmin = std::min(dx > 0. ? dx : 1e30, dy > 0. ? dy : 1e30);
		if (dmin > 0. && dmin < 1e29)
			nb = std::max(nb, (int)std::ceil(border_margin / dmin) + 1);
	}
	long nkeepc = 0;
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++) {
			const bool on_border =
				(border[0] && i <= nb - 1) ||            // left   (x-)
				(border[1] && i >= nx - nb) ||           // right  (x+)
				(border[2] && j <= nb - 1) ||            // bottom (y-)
				(border[3] && j >= ny - nb);             // top    (y+)
			IGNORE[i * ny + j] = on_border ? 0 : 1;
			if (on_border) nkeepc++;
		}
	// Depth of the stored strip in METRES is what a linked solver has to
	// compare against its own stencil reach; nodes alone do not say it.
	cfdwm_mpi::out() << "lsgrid: storing borders"
		<< (border[0] ? " left" : "") << (border[1] ? " right" : "")
		<< (border[2] ? " bottom" : "") << (border[3] ? " top" : "")
		<< " -- " << nb << " nodes + halo (margin "
		<< (border_margin > 0. ? border_margin : 0.) << " m), reaching "
		<< (nb - 1 + 2) * dx << " m into the domain in x and "
		<< (nb - 1 + 2) * dy << " m in y; " << nkeepc << " of "
		<< (long)nx * ny << " nodes before the stencil halo." << std::endl;
}
