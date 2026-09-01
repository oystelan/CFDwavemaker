#ifndef SEABED_H
#define SEABED_H

#include <string>
#include <vector>

/* Variable bathymetry from a triangulated surface file.

   The seabed is given as a geometric surface z_bed(x, y) (z up, still
   water level at z = 0); the water depth is

       h(x, y) = -(z_bed(x, y) + zshift)

   with h positive down. Supported formats (dispatch on file extension):
     .stl  - ascii or binary STL (auto-detected)
     .obj  - Wavefront OBJ (v/f records; polygons are fan-triangulated)

   Probing:
     - depth is evaluated on the facet planes (exact piecewise-linear
       interpolation of the mesh);
     - the gradient (dh/dx, dh/dy) comes from angle-weighted vertex
       normals interpolated barycentrically (Phong-style), so it is
       continuous across facet edges - important because the sloping-bed
       correction is linear in grad h, and facet-normal jumps would
       print straight into the kinematics;
     - where several triangles cover the same (x, y) (e.g. a closed
       solid exported from CAD), the highest surface is used;
     - queries outside the triangulated footprint clamp to the nearest
       point of the mesh (one-time warning) and return zero gradient,
       so a bathymetry patch smaller than the wave domain degrades
       gracefully to constant depth beyond its edge;
     - a minimum-depth limit clamps shallow regions: where the mesh is
       shallower than min_depth the reported depth is min_depth with
       zero gradient (a flat shelf). Potential-flow wave theory is
       invalid (and numerically explosive, 1/sinh(kh)) as h -> 0, so
       shoreline-crossing bathymetry must not be probed raw. The limit
       is a modelling parameter; set it from the sea state (waves break
       around h ~ 1.3*Hs).

   A 2D uniform bin index over the triangle xy-projections makes point
   queries O(1) after load. Triangles that are (near-)vertical in the
   xy-projection (side walls of a solid) are skipped. */
class Seabed {
public:
	Seabed() : active_(false), zshift_(0.), min_depth_(0.1), nbx_(0),
		nby_(0), bx0_(0.), by0_(0.), bdx_(1.), bdy_(1.),
		warned_outside_(false), warned_shallow_(false) {}

	// Load a surface file; returns false (with a message on stderr) on
	// failure. zshift moves the surface vertically before the depth
	// convention above is applied; xshift/yshift are subtracted from
	// the horizontal coordinates (x_domain = x_file - xshift), so
	// meshes with large coordinate offsets (e.g. UTM) can be
	// repositioned onto the wave domain. WaveForge auto-centers such
	// meshes and writes the shifts to the [seabed] section.
	bool load(const std::string& path, double zshift = 0.,
		double xshift = 0., double yshift = 0.);

	// kinematics gate: true when the surface is loaded AND the wave
	// kinematics should use it (slope corrections, local dispersion)
	bool active() const { return active_; }
	// surface present, independent of the kinematics gate. Geometry-only
	// mode ([seabed] kinematics 0) keeps loaded() true while active()
	// is false: the surface is still served to linked solvers through
	// wave_Seabed*() but the kinematics assume constant depth.
	bool loaded() const { return loaded_; }
	void set_kinematics(bool on) { active_ = loaded_ && on; }

	// Free the mesh once init-time sampling is done (geometry-only mode:
	// the kinematics never touch the seabed at runtime). Refused while
	// the kinematics use the surface (active()); any query after release
	// is a programming error and aborts loudly.
	void release();
	bool released() const { return released_; }

	// Minimum-depth clamp [m]: depth() never returns less than this and
	// the gradient is zeroed in clamped regions. Default 0.1 m.
	void set_min_depth(double hmin) { min_depth_ = hmin; }
	double get_min_depth() const { return min_depth_; }

	// water depth [m], positive down
	double depth(double x, double y) const;

	// unclamped water depth [m]: no minimum-depth limit, so land comes
	// out negative (bed above the surface z = -zshift plane). For
	// solvers importing the true bathymetry, NOT for wave kinematics.
	double depth_raw(double x, double y) const;

	// horizontal depth gradient (dh/dx, dh/dy)
	void gradient(double x, double y, double& hx, double& hy) const;

	// unit surface normal (pointing up, out of the bed)
	void normal(double x, double y, double n[3]) const;

	// depth range over the mesh vertices (after zshift)
	double min_depth() const { return hmin_; }
	double max_depth() const { return hmax_; }

	// xy bounding box of the triangulated footprint (after the shifts).
	// A solver must keep its whole domain inside this box AND inside the
	// footprint itself (the mesh need not fill its own bbox) -- queries
	// beyond the edge are silently clamped to the nearest mesh point,
	// which shows up as a flat, zero-slope skirt rather than an error.
	void bounds(double& xmin, double& xmax,
		double& ymin, double& ymax) const;

	// is (x, y) covered by an actual triangle?  No clamping, no warning:
	// the exact test a solver should use to validate its domain up front.
	bool contains(double x, double y) const;

	// how many queries have been served by out-of-footprint clamping so
	// far (see the warning in locate()).  Non-zero means part of the
	// sampled region is not really described by the bathymetry.
	long outside_count() const { return noutside_; }
	void reset_outside_count() const { noutside_ = 0; }

	size_t ntriangles() const { return tri_.size() / 3; }
	size_t nvertices() const { return vx_.size(); }

private:
	// parse + weld + finalize on the calling rank (the serial path; under
	// MPI only rank 0 runs it and the result is broadcast)
	bool load_local(const std::string& path, double xshift, double yshift);
	bool load_stl(const std::string& path);
	bool load_stl_ascii(const std::string& path);
	bool load_stl_binary(const std::string& path);
	bool load_obj(const std::string& path);
	// weld duplicated vertices, drop degenerate/vertical triangles,
	// orient windings upward, build vertex normals and the bin index
	bool finalize(const std::string& path);
	void build_bins();

	// locate the triangle covering (x, y); returns -1 if outside the
	// footprint. On success wgt holds the barycentric weights.
	int find_triangle(double x, double y, double wgt[3]) const;
	// find_triangle with out-of-footprint fallback to the nearest mesh
	// point, memoized on the last query point (the wave kinematics probe
	// the same (x, y) once per layer and order term). Sets outside.
	int locate(double x, double y, double wgt[3], bool& outside) const;
	// nearest point of the mesh to (x, y) in the xy-plane (for
	// out-of-footprint clamping); returns triangle index, fills wgt.
	int nearest_triangle(double x, double y, double wgt[3]) const;
	// interpolated upward unit normal for triangle t at weights wgt
	void interp_normal(int t, const double wgt[3], double n[3]) const;

	bool active_;
	bool loaded_ = false;
	bool released_ = false;
	double zshift_;
	double min_depth_;
	std::vector<double> vx_, vy_, vz_;    // welded vertices
	std::vector<int> tri_;                // 3 vertex indices per triangle
	std::vector<double> nx_, ny_, nz_;    // vertex normals (unit, up)
	double hmin_, hmax_;

	// uniform bin index over the xy bounding box
	int nbx_, nby_;
	double bx0_, by0_, bdx_, bdy_;
	std::vector<std::vector<int> > bins_;

	mutable bool warned_outside_;
	mutable bool warned_shallow_;
	mutable long noutside_ = 0;   // out-of-footprint queries served by clamping
	// unique id per successful load(); keys the thread-local locate memo
	// so a reloaded or reallocated Seabed can never serve stale hits
	int load_id_ = 0;
};

#endif
