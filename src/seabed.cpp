#include "seabed.h"
#include "mpi_utils.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <unordered_map>

namespace {

std::string lower_ext(const std::string& path) {
	size_t dot = path.find_last_of('.');
	if (dot == std::string::npos) return "";
	std::string ext = path.substr(dot + 1);
	for (size_t i = 0; i < ext.size(); i++)
		ext[i] = (char)std::tolower((unsigned char)ext[i]);
	return ext;
}

// squared distance from p to segment a-b in 2D; wab returns the
// parameter of the closest point along the segment (0 at a, 1 at b)
double seg_dist2(double px, double py, double ax, double ay,
                 double bx, double by, double& wab) {
	double ex = bx - ax, ey = by - ay;
	double L2 = ex * ex + ey * ey;
	double t = 0.;
	if (L2 > 0.)
		t = std::max(0., std::min(1., ((px - ax) * ex + (py - ay) * ey) / L2));
	double dx = px - (ax + t * ex), dy = py - (ay + t * ey);
	wab = t;
	return dx * dx + dy * dy;
}

} // namespace


bool Seabed::load(const std::string& path, double zshift,
	double xshift, double yshift) {
	active_ = false;
	zshift_ = zshift;
	vx_.clear(); vy_.clear(); vz_.clear();
	tri_.clear(); nx_.clear(); ny_.clear(); nz_.clear();
	bins_.clear();
	warned_outside_ = false;
	noutside_ = 0;

#if defined(MPI_enable)
	/* Parse and weld on rank 0 only, then broadcast the compact final
	   mesh. N ranks concurrently parsing a large ascii OBJ (transient
	   parse + weld buffers run ~3x the final arrays) is what used to
	   surge node memory at init. The broadcast arrays answer queries
	   exactly like a local load; the bin index is rebuilt locally
	   (cheap O(ntri)). */
	if (cfdwm_mpi::active() && cfdwm_mpi::size() > 1) {
		int root_ok = 0;
		if (cfdwm_mpi::rank() == 0)
			root_ok = load_local(path, xshift, yshift) ? 1 : 0;
		MPI_Bcast(&root_ok, 1, MPI_INT, 0, cfdwm_mpi::comm());
		if (!root_ok) return false;
		long n[2] = { (long)vx_.size(), (long)tri_.size() };
		MPI_Bcast(n, 2, MPI_LONG, 0, cfdwm_mpi::comm());
		if (cfdwm_mpi::rank() != 0) {
			vx_.resize(n[0]); vy_.resize(n[0]); vz_.resize(n[0]);
			nx_.resize(n[0]); ny_.resize(n[0]); nz_.resize(n[0]);
			tri_.resize(n[1]);
		}
		MPI_Bcast(vx_.data(), n[0], MPI_DOUBLE, 0, cfdwm_mpi::comm());
		MPI_Bcast(vy_.data(), n[0], MPI_DOUBLE, 0, cfdwm_mpi::comm());
		MPI_Bcast(vz_.data(), n[0], MPI_DOUBLE, 0, cfdwm_mpi::comm());
		MPI_Bcast(nx_.data(), n[0], MPI_DOUBLE, 0, cfdwm_mpi::comm());
		MPI_Bcast(ny_.data(), n[0], MPI_DOUBLE, 0, cfdwm_mpi::comm());
		MPI_Bcast(nz_.data(), n[0], MPI_DOUBLE, 0, cfdwm_mpi::comm());
		MPI_Bcast(tri_.data(), n[1], MPI_INT, 0, cfdwm_mpi::comm());
		double hh[2] = { hmin_, hmax_ };
		MPI_Bcast(hh, 2, MPI_DOUBLE, 0, cfdwm_mpi::comm());
		hmin_ = hh[0]; hmax_ = hh[1];
		if (cfdwm_mpi::rank() != 0)
			build_bins();
	}
	else if (!load_local(path, xshift, yshift))
		return false;
#else
	if (!load_local(path, xshift, yshift))
		return false;
#endif
	// unique per successful load (single-threaded init): invalidates any
	// thread-local locate memo pointing at a previous mesh
	static int load_seq = 0;
	load_id_ = ++load_seq;
	released_ = false;
	loaded_ = true;
	active_ = true;
	return true;
}

bool Seabed::load_local(const std::string& path, double xshift, double yshift) {
	std::string ext = lower_ext(path);
	bool ok;
	if (ext == "stl")
		ok = load_stl(path);
	else if (ext == "obj")
		ok = load_obj(path);
	else {
		fprintf(stderr, "Seabed: unsupported file format \".%s\" (%s). "
			"Supported: .stl, .obj\n", ext.c_str(), path.c_str());
		return false;
	}
	if (!ok) return false;
	if (xshift != 0. || yshift != 0.) {
		for (size_t i = 0; i < vx_.size(); i++) {
			vx_[i] -= xshift;
			vy_[i] -= yshift;
		}
	}
	return finalize(path);
}

/* Free the mesh once init-time sampling is done (geometry-only mode: the
   kinematics never touch the seabed at runtime). Refused while the
   kinematics use the surface; queries after release abort loudly. */
void Seabed::release() {
	if (!loaded_ || released_) return;
	if (active_) {
		fprintf(stderr, "Seabed: release() refused -- the wave kinematics "
			"use the surface ([seabed] kinematics 1).\n");
		return;
	}
	vx_.clear(); vx_.shrink_to_fit();
	vy_.clear(); vy_.shrink_to_fit();
	vz_.clear(); vz_.shrink_to_fit();
	nx_.clear(); nx_.shrink_to_fit();
	ny_.clear(); ny_.shrink_to_fit();
	nz_.clear(); nz_.shrink_to_fit();
	tri_.clear(); tri_.shrink_to_fit();
	bins_.clear(); bins_.shrink_to_fit();
	released_ = true;
}


// ---------------------------------------------------------------- STL

bool Seabed::load_stl(const std::string& path) {
	// A binary STL is exactly 84 + 50*ntri bytes; use that as the
	// primary discriminator (ascii files starting with "solid" exist,
	// and so do binary files whose 80-byte header starts with "solid").
	std::ifstream fh(path.c_str(), std::ios::binary | std::ios::ate);
	if (!fh) {
		fprintf(stderr, "Seabed: cannot open %s\n", path.c_str());
		return false;
	}
	std::streamoff size = fh.tellg();
	if (size >= 84) {
		fh.seekg(80);
		uint32_t ntri = 0;
		fh.read(reinterpret_cast<char*>(&ntri), 4);
		if (fh && size == (std::streamoff)84 + (std::streamoff)50 * ntri)
			return load_stl_binary(path);
	}
	return load_stl_ascii(path);
}

bool Seabed::load_stl_ascii(const std::string& path) {
	std::ifstream fh(path.c_str());
	if (!fh) {
		fprintf(stderr, "Seabed: cannot open %s\n", path.c_str());
		return false;
	}
	std::string tok;
	int nv = 0;
	double v[3];
	while (fh >> tok) {
		if (tok == "vertex") {
			if (!(fh >> v[0] >> v[1] >> v[2])) {
				fprintf(stderr, "Seabed: malformed vertex in %s\n",
					path.c_str());
				return false;
			}
			vx_.push_back(v[0]);
			vy_.push_back(v[1]);
			vz_.push_back(v[2]);
			if (++nv == 3) {
				int b = (int)vx_.size() - 3;
				tri_.push_back(b);
				tri_.push_back(b + 1);
				tri_.push_back(b + 2);
				nv = 0;
			}
		}
	}
	if (nv != 0) {
		fprintf(stderr, "Seabed: %s ends mid-facet\n", path.c_str());
		return false;
	}
	return true;
}

bool Seabed::load_stl_binary(const std::string& path) {
	std::ifstream fh(path.c_str(), std::ios::binary);
	if (!fh) {
		fprintf(stderr, "Seabed: cannot open %s\n", path.c_str());
		return false;
	}
	fh.seekg(80);
	uint32_t ntri = 0;
	fh.read(reinterpret_cast<char*>(&ntri), 4);
	vx_.reserve(3 * ntri); vy_.reserve(3 * ntri); vz_.reserve(3 * ntri);
	tri_.reserve(3 * ntri);
	for (uint32_t t = 0; t < ntri; t++) {
		float rec[12];
		fh.read(reinterpret_cast<char*>(rec), 48);
		fh.seekg(2, std::ios::cur); // attribute byte count
		if (!fh) {
			fprintf(stderr, "Seabed: truncated binary STL %s\n",
				path.c_str());
			return false;
		}
		int b = (int)vx_.size();
		for (int c = 0; c < 3; c++) { // rec[0..2] is the facet normal
			vx_.push_back(rec[3 + 3 * c]);
			vy_.push_back(rec[4 + 3 * c]);
			vz_.push_back(rec[5 + 3 * c]);
		}
		tri_.push_back(b);
		tri_.push_back(b + 1);
		tri_.push_back(b + 2);
	}
	return true;
}


// ---------------------------------------------------------------- OBJ

bool Seabed::load_obj(const std::string& path) {
	std::ifstream fh(path.c_str());
	if (!fh) {
		fprintf(stderr, "Seabed: cannot open %s\n", path.c_str());
		return false;
	}
	std::string line;
	while (std::getline(fh, line)) {
		std::istringstream ls(line);
		std::string key;
		if (!(ls >> key)) continue;
		if (key == "v") {
			double x, y, z;
			if (!(ls >> x >> y >> z)) {
				fprintf(stderr, "Seabed: malformed OBJ vertex in %s\n",
					path.c_str());
				return false;
			}
			vx_.push_back(x); vy_.push_back(y); vz_.push_back(z);
		} else if (key == "f") {
			std::vector<int> idx;
			std::string ref;
			while (ls >> ref) {
				// "i", "i/j", "i//k", "i/j/k"; 1-based, negative = from end
				size_t slash = ref.find('/');
				long i = std::atol(ref.substr(0, slash).c_str());
				if (i < 0) i += (long)vx_.size() + 1;
				if (i < 1 || i >(long)vx_.size()) {
					fprintf(stderr, "Seabed: OBJ face index out of range "
						"in %s\n", path.c_str());
					return false;
				}
				idx.push_back((int)(i - 1));
			}
			if (idx.size() < 3) {
				fprintf(stderr, "Seabed: OBJ face with <3 vertices in %s\n",
					path.c_str());
				return false;
			}
			for (size_t c = 1; c + 1 < idx.size(); c++) { // fan
				tri_.push_back(idx[0]);
				tri_.push_back(idx[c]);
				tri_.push_back(idx[c + 1]);
			}
		}
	}
	return true;
}


// ------------------------------------------------- weld, normals, bins

bool Seabed::finalize(const std::string& path) {
	size_t nraw = vx_.size();
	if (nraw < 3 || tri_.size() < 3) {
		fprintf(stderr, "Seabed: no triangles found in %s\n", path.c_str());
		return false;
	}

	// bounding box (raw) for the weld tolerance
	double xmin = vx_[0], xmax = vx_[0], ymin = vy_[0], ymax = vy_[0];
	double zmin = vz_[0], zmax = vz_[0];
	for (size_t i = 1; i < nraw; i++) {
		xmin = std::min(xmin, vx_[i]); xmax = std::max(xmax, vx_[i]);
		ymin = std::min(ymin, vy_[i]); ymax = std::max(ymax, vy_[i]);
		zmin = std::min(zmin, vz_[i]); zmax = std::max(zmax, vz_[i]);
	}
	double scale = std::max(std::max(xmax - xmin, ymax - ymin),
		std::max(zmax - zmin, 1e-30));
	double qtol = 1e-7 * scale;

	// weld duplicated vertices (STL repeats them per facet) so vertex
	// normals can be shared across facets
	std::unordered_map<std::string, int> vmap;
	std::vector<int> remap(nraw);
	std::vector<double> wx, wy, wz;
	wx.reserve(nraw / 4); wy.reserve(nraw / 4); wz.reserve(nraw / 4);
	for (size_t i = 0; i < nraw; i++) {
		int64_t q[3] = { (int64_t)std::llround(vx_[i] / qtol),
			(int64_t)std::llround(vy_[i] / qtol),
			(int64_t)std::llround(vz_[i] / qtol) };
		std::string key(reinterpret_cast<char*>(q), sizeof(q));
		std::unordered_map<std::string, int>::iterator it = vmap.find(key);
		if (it == vmap.end()) {
			int id = (int)wx.size();
			vmap.insert(std::make_pair(key, id));
			wx.push_back(vx_[i]); wy.push_back(vy_[i]); wz.push_back(vz_[i]);
			remap[i] = id;
		} else {
			remap[i] = it->second;
		}
	}
	vx_.swap(wx); vy_.swap(wy); vz_.swap(wz);

	// remap triangles; orient upward; drop triangles that are
	// degenerate or (near-)vertical in the xy-projection (side walls)
	std::vector<int> tris;
	tris.reserve(tri_.size());
	size_t nskipped = 0;
	for (size_t t = 0; t + 2 < tri_.size(); t += 3) {
		int a = remap[tri_[t]], b = remap[tri_[t + 1]], c = remap[tri_[t + 2]];
		if (a == b || b == c || a == c) { nskipped++; continue; }
		double e1x = vx_[b] - vx_[a], e1y = vy_[b] - vy_[a],
			e1z = vz_[b] - vz_[a];
		double e2x = vx_[c] - vx_[a], e2y = vy_[c] - vy_[a],
			e2z = vz_[c] - vz_[a];
		double nx = e1y * e2z - e1z * e2y;
		double ny = e1z * e2x - e1x * e2z;
		double nz = e1x * e2y - e1y * e2x; // = 2 * signed xy-area
		double nn = std::sqrt(nx * nx + ny * ny + nz * nz);
		if (nn <= 0. || std::fabs(nz) < 1e-8 * nn) { nskipped++; continue; }
		if (nz < 0.) std::swap(b, c); // wind upward
		tris.push_back(a); tris.push_back(b); tris.push_back(c);
	}
	tri_.swap(tris);
	if (tri_.empty()) {
		fprintf(stderr, "Seabed: %s has no horizontally-projectable "
			"triangles\n", path.c_str());
		return false;
	}
	if (nskipped)
		printf("Seabed: skipped %zu degenerate/vertical facets in %s\n",
			nskipped, path.c_str());

	// angle-weighted vertex normals (all triangles now wind upward)
	size_t nv = vx_.size();
	nx_.assign(nv, 0.); ny_.assign(nv, 0.); nz_.assign(nv, 0.);
	for (size_t t = 0; t < tri_.size(); t += 3) {
		int id[3] = { tri_[t], tri_[t + 1], tri_[t + 2] };
		double fnx, fny, fnz;
		{
			double e1x = vx_[id[1]] - vx_[id[0]],
				e1y = vy_[id[1]] - vy_[id[0]],
				e1z = vz_[id[1]] - vz_[id[0]];
			double e2x = vx_[id[2]] - vx_[id[0]],
				e2y = vy_[id[2]] - vy_[id[0]],
				e2z = vz_[id[2]] - vz_[id[0]];
			fnx = e1y * e2z - e1z * e2y;
			fny = e1z * e2x - e1x * e2z;
			fnz = e1x * e2y - e1y * e2x;
			double nn = std::sqrt(fnx * fnx + fny * fny + fnz * fnz);
			fnx /= nn; fny /= nn; fnz /= nn;
		}
		for (int c = 0; c < 3; c++) {
			int i0 = id[c], i1 = id[(c + 1) % 3], i2 = id[(c + 2) % 3];
			double ax = vx_[i1] - vx_[i0], ay = vy_[i1] - vy_[i0],
				az = vz_[i1] - vz_[i0];
			double bx = vx_[i2] - vx_[i0], by = vy_[i2] - vy_[i0],
				bz = vz_[i2] - vz_[i0];
			double la = std::sqrt(ax * ax + ay * ay + az * az);
			double lb = std::sqrt(bx * bx + by * by + bz * bz);
			double ca = (ax * bx + ay * by + az * bz) / (la * lb);
			ca = std::max(-1., std::min(1., ca));
			double ang = std::acos(ca);
			nx_[i0] += ang * fnx;
			ny_[i0] += ang * fny;
			nz_[i0] += ang * fnz;
		}
	}
	for (size_t i = 0; i < nv; i++) {
		double nn = std::sqrt(nx_[i] * nx_[i] + ny_[i] * ny_[i]
			+ nz_[i] * nz_[i]);
		if (nn > 0.) { nx_[i] /= nn; ny_[i] /= nn; nz_[i] /= nn; }
		else { nx_[i] = ny_[i] = 0.; nz_[i] = 1.; }
	}

	hmin_ = hmax_ = -(vz_[tri_[0]] + zshift_);
	for (size_t t = 0; t < tri_.size(); t++) {
		double h = -(vz_[tri_[t]] + zshift_);
		hmin_ = std::min(hmin_, h);
		hmax_ = std::max(hmax_, h);
	}

	build_bins();
	return true;
}

void Seabed::build_bins() {
	size_t ntri = tri_.size() / 3;
	bx0_ = vx_[tri_[0]]; by0_ = vy_[tri_[0]];
	double x1 = bx0_, y1 = by0_;
	for (size_t t = 0; t < tri_.size(); t++) {
		bx0_ = std::min(bx0_, vx_[tri_[t]]);
		x1 = std::max(x1, vx_[tri_[t]]);
		by0_ = std::min(by0_, vy_[tri_[t]]);
		y1 = std::max(y1, vy_[tri_[t]]);
	}
	int nb = (int)std::ceil(std::sqrt((double)std::max(ntri, (size_t)1) / 2.));
	nbx_ = std::max(8, std::min(2048, nb));
	nby_ = nbx_;
	bdx_ = std::max((x1 - bx0_) / nbx_, 1e-30);
	bdy_ = std::max((y1 - by0_) / nby_, 1e-30);
	bins_.assign((size_t)nbx_ * nby_, std::vector<int>());
	for (size_t t = 0; t < ntri; t++) {
		int a = tri_[3 * t], b = tri_[3 * t + 1], c = tri_[3 * t + 2];
		double txmin = std::min(vx_[a], std::min(vx_[b], vx_[c]));
		double txmax = std::max(vx_[a], std::max(vx_[b], vx_[c]));
		double tymin = std::min(vy_[a], std::min(vy_[b], vy_[c]));
		double tymax = std::max(vy_[a], std::max(vy_[b], vy_[c]));
		int i0 = std::max(0, std::min(nbx_ - 1,
			(int)((txmin - bx0_) / bdx_)));
		int i1 = std::max(0, std::min(nbx_ - 1,
			(int)((txmax - bx0_) / bdx_)));
		int j0 = std::max(0, std::min(nby_ - 1,
			(int)((tymin - by0_) / bdy_)));
		int j1 = std::max(0, std::min(nby_ - 1,
			(int)((tymax - by0_) / bdy_)));
		for (int i = i0; i <= i1; i++)
			for (int j = j0; j <= j1; j++)
				bins_[(size_t)i * nby_ + j].push_back((int)t);
	}
}


// ----------------------------------------------------------- queries

int Seabed::find_triangle(double x, double y, double wgt[3]) const {
	int bi = (int)((x - bx0_) / bdx_), bj = (int)((y - by0_) / bdy_);
	if (bi < 0 || bi >= nbx_ || bj < 0 || bj >= nby_) return -1;
	const std::vector<int>& cand = bins_[(size_t)bi * nby_ + bj];
	int best = -1;
	double bestz = 0., bw[3];
	for (size_t q = 0; q < cand.size(); q++) {
		int t = cand[q];
		int a = tri_[3 * t], b = tri_[3 * t + 1], c = tri_[3 * t + 2];
		double d = (vx_[b] - vx_[a]) * (vy_[c] - vy_[a])
			- (vx_[c] - vx_[a]) * (vy_[b] - vy_[a]);
		// d > 0 by construction (upward winding, no degenerates)
		double w1 = ((x - vx_[a]) * (vy_[c] - vy_[a])
			- (vx_[c] - vx_[a]) * (y - vy_[a])) / d;
		double w2 = ((vx_[b] - vx_[a]) * (y - vy_[a])
			- (x - vx_[a]) * (vy_[b] - vy_[a])) / d;
		double w0 = 1. - w1 - w2;
		double eps = 1e-10;
		if (w0 < -eps || w1 < -eps || w2 < -eps) continue;
		double z = w0 * vz_[a] + w1 * vz_[b] + w2 * vz_[c];
		// several triangles may cover this (x, y) (closed solid):
		// the seabed is the highest surface
		if (best < 0 || z > bestz) {
			best = t; bestz = z;
			bw[0] = w0; bw[1] = w1; bw[2] = w2;
		}
	}
	if (best >= 0) { wgt[0] = bw[0]; wgt[1] = bw[1]; wgt[2] = bw[2]; }
	return best;
}

int Seabed::nearest_triangle(double x, double y, double wgt[3]) const {
	int bi = std::max(0, std::min(nbx_ - 1, (int)((x - bx0_) / bdx_)));
	int bj = std::max(0, std::min(nby_ - 1, (int)((y - by0_) / bdy_)));
	// distance from the query to the bin grid (zero inside): all mesh
	// geometry lies within the grid, so this offset adds in quadrature
	// to the lateral ring distance in the stopping bound below. Without
	// it, a far-outside query scans every ring within its full clamp
	// distance (cost growing quadratically with the overhang).
	double cx = std::min(std::max(x, bx0_), bx0_ + nbx_ * bdx_);
	double cy = std::min(std::max(y, by0_), by0_ + nby_ * bdy_);
	double dout2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
	int best = -1;
	double bestd2 = 0.;
	int maxr = std::max(nbx_, nby_);
	for (int r = 0; r <= maxr; r++) {
		if (best >= 0) {
			// stop once closer bins are exhausted: bins at ring r are
			// at least (r-1)*min(bdx,bdy) away laterally
			double ringmin = (r - 1.) * std::min(bdx_, bdy_);
			if (ringmin > 0. && dout2 + ringmin * ringmin > bestd2) break;
		}
		for (int i = bi - r; i <= bi + r; i++) {
			for (int j = bj - r; j <= bj + r; j++) {
				if (std::max(std::abs(i - bi), std::abs(j - bj)) != r)
					continue; // ring only
				if (i < 0 || i >= nbx_ || j < 0 || j >= nby_) continue;
				const std::vector<int>& cand = bins_[(size_t)i * nby_ + j];
				for (size_t q = 0; q < cand.size(); q++) {
					int t = cand[q];
					int a = tri_[3 * t], b = tri_[3 * t + 1],
						c = tri_[3 * t + 2];
					// distance to the triangle boundary in 2D
					const int va[3] = { a, b, c }, vb[3] = { b, c, a };
					for (int e = 0; e < 3; e++) {
						double s;
						double d2 = seg_dist2(x, y, vx_[va[e]], vy_[va[e]],
							vx_[vb[e]], vy_[vb[e]], s);
						if (best < 0 || d2 < bestd2) {
							best = t; bestd2 = d2;
							wgt[0] = wgt[1] = wgt[2] = 0.;
							int ia = (e == 0) ? 0 : (e == 1) ? 1 : 2;
							int ib = (ia + 1) % 3;
							wgt[ia] = 1. - s;
							wgt[ib] = s;
						}
					}
				}
			}
		}
	}
	return best;
}

void Seabed::interp_normal(int t, const double wgt[3], double n[3]) const {
	int a = tri_[3 * t], b = tri_[3 * t + 1], c = tri_[3 * t + 2];
	n[0] = wgt[0] * nx_[a] + wgt[1] * nx_[b] + wgt[2] * nx_[c];
	n[1] = wgt[0] * ny_[a] + wgt[1] * ny_[b] + wgt[2] * ny_[c];
	n[2] = wgt[0] * nz_[a] + wgt[1] * nz_[b] + wgt[2] * nz_[c];
	double nn = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	if (nn > 0.) { n[0] /= nn; n[1] /= nn; n[2] /= nn; }
	else { n[0] = n[1] = 0.; n[2] = 1.; }
}

namespace {
// Last-query memo for Seabed::locate. The wave kinematics probe depth
// and gradient at the same (x, y) once per layer and per order term
// (~6 probes per grid point-layer), so a single-entry memo removes
// nearly all repeat triangle searches. thread_local: OpenMP threads
// each keep their own entry; keyed on load_id so a reloaded (or
// reallocated) Seabed can never serve a stale hit.
struct LocateMemo {
	int load_id;
	double x, y;
	int t;
	bool outside;
	double w[3];
};
thread_local LocateMemo locate_memo = { 0, 0., 0., -1, false,
	{ 0., 0., 0. } };
} // namespace

int Seabed::locate(double x, double y, double wgt[3], bool& outside) const {
	if (released_) {
		fprintf(stderr, "[rank %d] Seabed ERROR: query after release() -- "
			"the mesh was freed after init-time sampling. Remove the "
			"wave_SeabedRelease() call if runtime queries are needed.\n",
			cfdwm_mpi::rank());
		cfdwm_mpi::abort(1);
	}

	LocateMemo& m = locate_memo;
	if (m.load_id == load_id_ && m.x == x && m.y == y) {
		wgt[0] = m.w[0]; wgt[1] = m.w[1]; wgt[2] = m.w[2];
		outside = m.outside;
		return m.t;
	}
	int t = find_triangle(x, y, wgt);
	outside = (t < 0);
	if (outside) {
		noutside_++;
		if (!warned_outside_) {
			warned_outside_ = true;
			// cfdwm_mpi::err(), not printf: a bare printf escapes both
			// basilisk's per-rank stdout redirection and the master-only
			// logging, so it used to print from every rank at once.
			cfdwm_mpi::err()
				<< "Seabed WARNING: query point (" << x << ", " << y
				<< ") is outside the bathymetry footprint; clamping to the"
				" nearest mesh point (constant depth, zero slope beyond"
				" the edge). Further occurrences are counted, not printed"
				" -- see wave_SeabedOutsideCount(). A domain reaching the"
				" edge of the bathymetry is a common cause of boundary"
				" blow-up." << std::endl;
		}
		t = nearest_triangle(x, y, wgt);
	}
	m.load_id = load_id_;
	m.x = x; m.y = y;
	m.t = t;
	m.outside = outside;
	m.w[0] = wgt[0]; m.w[1] = wgt[1]; m.w[2] = wgt[2];
	return t;
}

void Seabed::bounds(double& xmin, double& xmax,
	double& ymin, double& ymax) const {
	// the bin index is built over the xy bounding box of the triangles,
	// so its extent IS the footprint bbox (already shifted)
	if (!loaded_ || bins_.empty()) {
		xmin = xmax = ymin = ymax = 0.;
		return;
	}
	xmin = bx0_;
	xmax = bx0_ + nbx_ * bdx_;
	ymin = by0_;
	ymax = by0_ + nby_ * bdy_;
}

bool Seabed::contains(double x, double y) const {
	if (!loaded_) return false;
	double w[3];
	return find_triangle(x, y, w) >= 0;
}

double Seabed::depth_raw(double x, double y) const {
	double w[3];
	bool outside;
	int t = locate(x, y, w, outside);
	if (t < 0) return 0.; // cannot happen with a non-empty mesh
	int a = tri_[3 * t], b = tri_[3 * t + 1], c = tri_[3 * t + 2];
	double z = w[0] * vz_[a] + w[1] * vz_[b] + w[2] * vz_[c];
	return -(z + zshift_);
}

double Seabed::depth(double x, double y) const {
	double h = depth_raw(x, y);
	if (h < min_depth_) {
		if (!warned_shallow_) {
			warned_shallow_ = true;
			printf("Seabed WARNING: bathymetry shallower than the "
				"minimum-depth limit %g m near (%g, %g); depth is "
				"clamped there (flat shelf, zero slope). Potential "
				"theory is invalid as h -> 0. Further occurrences "
				"are silent.\n", min_depth_, x, y);
		}
		return min_depth_;
	}
	return h;
}

void Seabed::gradient(double x, double y, double& hx, double& hy) const {
	double w[3];
	bool outside;
	int t = locate(x, y, w, outside);
	if (t < 0 || outside) { hx = hy = 0.; return; } // outside: constant depth
	int a = tri_[3 * t], b = tri_[3 * t + 1], c = tri_[3 * t + 2];
	double z = w[0] * vz_[a] + w[1] * vz_[b] + w[2] * vz_[c];
	if (-(z + zshift_) < min_depth_) { hx = hy = 0.; return; } // clamped
	double n[3];
	interp_normal(t, w, n);
	double nz = std::max(n[2], 1e-8);
	// h = -(z + zshift)  =>  grad h = -grad z = (nx, ny)/nz
	hx = n[0] / nz;
	hy = n[1] / nz;
}

void Seabed::normal(double x, double y, double n[3]) const {
	double w[3];
	bool outside;
	int t = locate(x, y, w, outside);
	if (t < 0 || outside) { n[0] = n[1] = 0.; n[2] = 1.; return; }
	int a = tri_[3 * t], b = tri_[3 * t + 1], c = tri_[3 * t + 2];
	double z = w[0] * vz_[a] + w[1] * vz_[b] + w[2] * vz_[c];
	if (-(z + zshift_) < min_depth_) { // clamped: flat shelf
		n[0] = n[1] = 0.; n[2] = 1.;
		return;
	}
	interp_normal(t, w, n);
}
