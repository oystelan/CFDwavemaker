#ifndef CFDWM_MPI_UTILS_H
#define CFDWM_MPI_UTILS_H

// Optional MPI support for distributing the lsGridSpline second-order kinematics
// computation across ranks. Guarded by -DMPI_enable=1; when the macro is absent
// the helpers collapse to a trivial serial implementation (rank 0 of 1), so the
// rest of the code is written once and compiles/runs unchanged without MPI.
//
// Contract: the host CFD owns MPI (it calls MPI_Init/Finalize). This library only
// *uses* the existing context; it never initialises or finalises MPI. All ranks in
// the communicator must reach the collective grid update together (which they do,
// because sgrids.update() is only driven by the collective wave_force_update()).

#include <algorithm>
#include <vector>

#if defined(MPI_enable)
#include <mpi.h>

namespace cfdwm_mpi {

	// Communicator used for the grid decomposition. Defaults to MPI_COMM_WORLD;
	// a host on a sub-communicator can override it via set_comm().
	inline MPI_Comm& comm() { static MPI_Comm c = MPI_COMM_WORLD; return c; }
	inline void set_comm(MPI_Comm c) { comm() = c; }

	inline bool active() { int f = 0; MPI_Initialized(&f); return f != 0; }
	inline int rank() { if (!active()) return 0; int r; MPI_Comm_rank(comm(), &r); return r; }
	inline int size() { if (!active()) return 1; int s; MPI_Comm_size(comm(), &s); return s; }

	// Even 1-D block decomposition of n items across the ranks: sets [i0,i1) for a
	// given rank r (larger remainder chunks go to the lowest ranks).
	inline void range_for(int n, int r, int s, int& i0, int& i1) {
		int chunk = n / s, rem = n % s;
		i0 = r * chunk + std::min(r, rem);
		i1 = i0 + chunk + (r < rem ? 1 : 0);
	}
	// [i0,i1) owned by *this* rank.
	inline void range(int n, int& i0, int& i1) { range_for(n, rank(), size(), i0, i1); }

	// Work-balanced contiguous partition of n columns into size() blocks, given a
	// per-column weight w[i] (>= 0, e.g. the number of active cells in the column).
	// Fills bounds[0..s] with bounds[0]=0, bounds[s]=n; deterministic on all ranks
	// (they all hold the same weights, so no communication is needed). Rank r owns
	// columns [bounds[r], bounds[r+1]).
	inline void partition_weighted(const std::vector<long>& w, int n, std::vector<int>& bounds) {
		int s = size();
		bounds.assign(s + 1, 0);
		bounds[s] = n;
		if (s <= 1) return;
		long total = 0; for (int i = 0; i < n; i++) total += w[i];
		if (total <= 0) { // degenerate (all-ignored): fall back to an even split
			for (int r = 0; r < s; r++) { int a, b; range_for(n, r, s, a, b); bounds[r] = a; }
			bounds[s] = n; return;
		}
		long acc = 0; int r = 1;
		for (int i = 0; i < n && r < s; i++) {
			acc += w[i];
			if (acc * s >= (long)r * total) { bounds[r] = i + 1; r++; }
		}
		for (; r < s; r++) bounds[r] = n;
		for (int rr = 1; rr < s; rr++) if (bounds[rr] < bounds[rr - 1]) bounds[rr] = bounds[rr - 1];
		bounds[s] = n;
	}

	// In-place all-gather of a slice laid out as [n][per_col] doubles, decomposed by
	// `bounds` (rank r owns columns [bounds[r],bounds[r+1])). Each rank has already
	// written its own block into `base`; afterwards every rank holds the full slice.
	inline void allgather_columns(double* base, const std::vector<int>& bounds, long per_col) {
		int s = size();
		if (s <= 1 || !active()) return;
		std::vector<int> counts(s), displs(s);
		for (int r = 0; r < s; r++) {
			displs[r] = (int)((long)bounds[r] * per_col);
			counts[r] = (int)((long)(bounds[r + 1] - bounds[r]) * per_col);
		}
		MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
		               base, counts.data(), displs.data(), MPI_DOUBLE, comm());
	}

} // namespace cfdwm_mpi

#else // ---- serial fallback (no MPI) ----------------------------------------

namespace cfdwm_mpi {
	inline int  rank() { return 0; }
	inline int  size() { return 1; }
	inline bool active() { return false; }
	inline void range(int n, int& i0, int& i1) { i0 = 0; i1 = n; }
	inline void partition_weighted(const std::vector<long>&, int n, std::vector<int>& bounds) {
		bounds.assign(2, 0); bounds[1] = n;
	}
	inline void allgather_columns(double*, const std::vector<int>&, long) {}
}

#endif

#endif // CFDWM_MPI_UTILS_H
