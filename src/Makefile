MAKE_DIR = $(PWD)

.DEFAULT_GOAL := default
.PHONY: default

clean:
	@$(MAKE) -f make_openmp.mk clean
	@$(MAKE) -f make_swd.mk clean
	@$(MAKE) -f make_all.mk clean

# build CFDwavemaker_swd_openmp (base)
basic:
	@$(MAKE) -f make_openmp.mk

# build CFDwavemaker_swd_openmp (base + swd)
default:
	@$(MAKE) -f make_swd.mk

# build entire CFDwavemaker_all (base + swd + vtk)
all:
	@$(MAKE) -f make_all.mk
