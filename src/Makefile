MAKE_DIR = $(PWD)

.DEFAULT_GOAL := default
.PHONY: default

clean:
	@$(MAKE) -f make_openmp.mk clean
	@$(MAKE) -f make_swd.mk clean
	@$(MAKE) -f make_all.mk clean

# build CFDwavemaker_openmp (base, no swd)
basic:
	@$(MAKE) -f make_openmp.mk

# build CFDwavemaker_mpi (base + MPI, distributes lsgrid kinematics across ranks)
mpi:
	@$(MAKE) -f make_mpi.mk clean
	@$(MAKE) -f make_mpi.mk

# build CFDwavemaker_swd (base + swd, openmp) -- the default
default:
	@$(MAKE) -f make_swd.mk

# build entire CFDwavemaker_all (base + swd + vtk)
all:
	@$(MAKE) -f make_all.mk
