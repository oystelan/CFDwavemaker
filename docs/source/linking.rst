Linking to your software
========================

Feedback has thought me that what people often find difficult is to link the library to their software of choice.
To hopefully help with this, a few examples of linking are provided below.

General
-------

If you have built the library (or downloaded the binaries), there are only 2 files you really need:

- `CFDwavemaker.h`: header file, tells you which functions are available
- `libCFDwavemaker_swd_openmp.a(or .so)`: the actually library

The file with the extension .a is the static link library, while the .so file is the dynamic link library.
`CFDwavemaker.h` is found in the ./src/ directory, while the libraries are found in the ./build/linux64 library.

.. note::

	Depending on which version of the library you have built, name of the .so/.a file will differ. To elaborate:

	- libCFDwavemaker_swd_openmp is the default build (built with swd extension, using openmp parallelization)
	- libCFDwavemaker_openmp is the basic build (built without swd extension (core library only), using openmp parallelization)
	- libCFDwavemaker_all_openmp is the complete build (built with swd and vtk extensions, using openmp parallelization)


Linking with C/C++
------------------

Linking up the CFDwavemaker to any C/C++ based code is fairly straight forward, and is done in exactly the same way as any other library you might want to add to your program.
Essentially you need to do a couple of things:

- make sure to add "#include CFDwavemaker.h" somewhere in the top of your program
- be sure to link up either the dynamic or the static library duing linking. If the libCFDwavemaker_openmp.so file is located in the same directory as your program, then adding a "-lCFDwavemaker_openmp" in your linking statement of your make file should do the job.
- if your compiling a native C program, then make sure to add "-lstdc++" when linking. 

A simple example is shown below (basilisk example)


Linking up with Python
----------------------

Occationally it may be useful generate kinematics data for other reasons than pushing it into a CFD solver. In this case python is a good and easy alternative, using the ctypes package. An example is given below.

.. code:: python

	from ctypes import*
	import matplotlib.pylab as plt
	import numpy as np
	import os

	# give location of dll
	#mydll = cdll.comflow
	temp = os.path.abspath(__file__)
	temp = os.path.realpath(temp)
	temp = os.path.dirname(temp)

	# we load the shared library 
	mydll = cdll.LoadLibrary("libCFDwavemaker_openmp.so")     

	# and define the functions which CFDwavemaker.h tells us are in the library
	def init_dll(mydll):
		aa = mydll.Init
		aa.restype = c_int
		aa.argtypes = [POINTER(c_double),POINTER(c_double)]
		return aa(c_double(0),c_double(0))

	def clean_up(mydll):
		aa = mydll.Cleanup
		aa.restype = c_int
		return aa()

	def wave_VeloX(mydll,t,x,y,z):
		aa = mydll.wave_VeloX
		aa.restype = c_double
		aa.argtypes = [c_double,c_double,c_double,c_double]
		return aa(c_double(x),c_double(y),c_double(z),c_double(t))

	def wave_VeloY(mydll,t,x,y,z):
		aa = mydll.wave_VeloY
		aa.restype = c_double
		aa.argtypes = [c_double,c_double,c_double,c_double]
		return aa(c_double(x),c_double(y),c_double(z),c_double(t))

	def wave_VeloZ(mydll,t,x,y,z):
		aa = mydll.wave_VeloZ
		aa.restype = c_double
		aa.argtypes = [c_double,c_double,c_double,c_double]
		return aa(c_double(x),c_double(y),c_double(z),c_double(t))

	def wave_SurfElev(mydll,t,x,y):
		aa = mydll.wave_SurfElev
		aa.restype = c_double
		aa.argtypes = [c_double,c_double,c_double]
		return aa(c_double(x),c_double(y),c_double(t))

	# All done. lets test:

	# we start by initializing (reads data from waveinput.dat) 
	init_dll(mydll)

	x = 4.8
	y = 0.23
	z = -0.05


	time = np.arange(0,20,0.05)

	elev = []
	u = []
	v = []
	w = []

	for t in time:
		elev.append(wave_SurfElev(mydll,t,x,y))
		u.append(wave_VeloX(mydll,t,x,y,z))
		v.append(wave_VeloY(mydll,t,x,y,z))
		w.append(wave_VeloZ(mydll,t,x,y,z))
		

	# now, lets view the resulting surface elevation
	plt.plot(time, elev, label="surface elevation")
	plt.legend()
	plt.xlabel('Time [sec]')
	plt.ylabel('Surface elevation [m]')
	plt.grid(True)
	plt.savefig("./result_eta.png")

	# and also the wave kinematics
	plt.clf()
	plt.plot(time,u, label='u')
	plt.plot(time,v, label='v')
	plt.plot(time,w, label='w')
	plt.xlabel('Time [sec]')
	plt.ylabel('particle velocity [m/s]')
	plt.legend()
	plt.grid(True)
	plt.savefig("./result_uvw.png")


	# all done. remember to clean up after us.
	clean_up(mydll)


Linking up with ComFLOW
-----------------------

`ComFLOW`_ is a Volume-of-fluids (VOF) Navier-Stokes solver for free-surface flow, and is an excellent CFD solver for modelling free surface waves. External libraries such as CFDwavemaker may be linked to the solver for providing kinematics using a predefined extern C function. These functions are available in CFDwavemaker, and therefore it is straight forward to use the library with ComFLOW once the shared library has been built. The instructions on how to do so is given below.

.. _`ComFLOW`: http://www.math.rug.nl/~veldman/comflow/comflow.html

1. Start by copying the CFDwavemaker.so library to some place ComFLOW can locate it. A good location is among the external library files located in the directory `$(COMFLOW_INSTALL_DIR)/lib/linux/`

2. In your `comflow.cfi` file (main control file for comflow simulation) specify that the external library should be used as show in the example xml code (extract from a `comflow.cfi` file) below.

3. Your done. Now you should be able to run ComFLOW using CFDwavemaker for initialization. Remember to provide `waveinput.dat` file in you comflow run folder when starting a simulation (one level up of the input_files folder). 

.. code:: xml

	...

	<!-- WAVES: Definition of incoming/initial wave field and current -->
      <waves start_with_still_water="true" initialize_fs="true" mean_depth="1.2">

      <!-- MODEL: Definition of wave model -->
           <model model="none"/>

      <!-- CURRENT: Current -->
           <current magnitude="0." angle="0."/>

      <!-- RAMPING: Ramping for smooth startup of a simulation -->
           <ramping ramptype="0" rampfs="1" ramp1="0." ramp2="0."/>
    </waves>
	<!-- COUPLING: Settings for coupling to a.o. moving objects (prescribed and 
      interactive motions), XMF mooring module, external solutions (e.g. other 
      ComFLOW simulations), ... -->
      <coupling>

      <!-- EXTERNAL_SOLUTION: Settings for coupling to a shared library -->
           <external_solution dllfile="CFDwavemaker_omp.so" initialize="true" ramp="false"/>

    </coupling>

    ...

If you wish to use CFDwavemaker for providing waves at the boundary, this is done by altering the `comflow.cfi` file. Reference is made to the `ComFLOW manual`_.

.. _`ComFLOW manual`: http://poseidon.housing.rug.nl/sphinx/

Linking up with Basilisk (or C)
-------------------------------

`Basilisk`_ is an open source library for the solution of partial differential equations on adaptive Cartesian meshes. The code is built in C (C99) with more than a few improvements to syntax (refered to as `Basilisk C`), making it efficient and fairly easy to use and understand. 
Linking CFDwavemaker to this library is straight forward, since it is built on the same programming language (almost). 

.. _`Basilisk`: http://www.basilisk.fr


This simple example makes use of the non-hydrostatic multilayer solver to propagate irregular waves in a domain which measures 1738x1738m. the water depth is set to 76.4m.

.. code:: C

	/**
	# Irregular wave case using periodic domain (multilayer solver)

	outputs .vts files which can be viewed in paraview
	The solution is obtained using the layered model and demonstrates its
	robustness and a degree of realism even for this complex case. 


	Made by: Oystein lande 2022  (modified version of basilisk/src/test/breaking.c)

	*/

	//#include "grid/multigrid.h"
	#include "layered/hydro.h"
	#include "layered/nh.h"
	#include "layered/remap.h"
	#include "layered/perfs.h"
	#include "view.h"


	#include "CFDwavemaker.h"

	/**
	The initial conditions is set in the waveinput.dat file and read through
	the CFDwavemaker lib. This particular example is a spatial periodic
	solution with a domain size L 1738m. We run the simulation for approximately
	5 mean wave periods.*/

	#define l_  1738.
	#define k_  (2.*pi)/l_
	#define h_  76.4
	#define g_  9.81
	#define T0  15.
	#define Tmax 5*T0

	/**
	The domain is periodic in $x$ and resolved using 128$^2$
	points and 20 layers. */

	int main()
	{
		size(l_);
		//omp_set_num_threads(1);
		origin(-L0 / 2., -L0 / 2.);
		periodic(right);
		periodic(top);
		N = 128;
		nl = 20;
		G = g_;
		//nu = 1. / RE;
		nu = 0;

		/* Initialize CFDwavemaker. waveinput.dat file is read when calling this function*/
		wave_Initialize();

		/**
		Some implicit damping is necessary to damp fast modes. This may be
		related to the slow/fast decoupling of the $\theta$-scheme mentioned
		by [Durran & Blossey, 2012](#durran2012). */

		//theta_H = 0.51;

		run();

		/*Release memory allocated to CFDwavemaker after simulation end*/
		wave_Cleanup();
	}

	/**
	The initial conditions for the free-surface and velocity are given by
	the third-order Stokes solution. */

	event init(i = 0)
	{

		/**
		We can use a larger CFL, in particular because we are not dealing
		with shallow-water/wetting/drying. */

		CFL = 0.8;

		/**
		The layer thicknesses follow a geometric progression, starting from
		a top layer with a thickness of 1/3 that of the regular
		distribution. */

		geometric_beta((1./3) * h_ / nl, true);


		// We set the seabed reference (zb), layer heights (h) and velocities (u.x u.y and w)
		foreach() {
			zb[] = -h_;
			double H = wave_SurfElev(x, y, 0) - zb[];
			double z = zb[];
			foreach_layer() {
				h[] = H * beta[point.l];
				z += h[] / 2.;
				u.x[] = wave_VeloX(x, y, z, 0);
				u.y[] = wave_VeloY(x, y, z, 0);
				w[] = wave_VeloZ(x, y, z, 0);
				z += h[] / 2.;
			}
		}
		boundary(all);
	}

	/**
	We add (an approximation of) horizontal viscosity. */

	event viscous_term(i++)
	horizontal_diffusion((scalar*) {u}, nu, dt);

	/**
	We log the evolution of the kinetic and potential energies.

	~~~gnuplot Evolution of the kinetic, potential and total energy
	set xlabel 't/T0'
	plot [0:6]'log' u 1:2 w l t 'kinetic', '' u 1:3 w l t 'potential', \
		'' u 1:(($2+$3)/2.) w l t 'total/2'
	~~~
	*/

	event logfile(i++; t <= Tmax)
	{
		double ke = 0., gpe = 0.;
		foreach(reduction(+:ke) reduction(+:gpe)) {
			foreach_layer() {
				double norm2 = sq(w[]);
				foreach_dimension()
					norm2 += sq(u.x[]);
				ke += norm2 * h[] * dv();
			}
			gpe += sq(eta[]) * dv();
		}
		fprintf(stderr, "%g %g %g\n", t / T0, ke / 2., g_ * gpe / 2.);
	}


	/**
	And generate the movie of the free surface, coloring with the horizontal
	particle velocity u.x */


	event movie(t += 0.1)
	{
		view(fov = 17.3106, quat = { 0.475152,0.161235,0.235565,0.832313 },
			tx = -0.0221727, ty = -0.0140227, width = 1200, height = 768);
		char s[80];
		sprintf(s, "t = %.2f T0", t / T0);
		draw_string(s, size = 80);
		//for (double x = -1; x <= 1; x++)
			//translate(x);
		squares("u.x", linear = true, z = "eta", min = -10., max = 12.);
		save("movie.mp4");
	}


The example may be compiles with the following `Makefile`

.. code:: bash

	KDT_LIBS=$(BASILISK)/kdt
	PPR_LIBS=$(BASILISK)/ppr
	GL_LIBS=$(BASILISK)/gl

	.PHONY: default clean

	default: irregular ;

	irregular: irregular.c
		CC99='gcc -std=c99' qcc irregular.c -o irregular -fopenmp -lCFDwavemaker_openmp -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran
		
	clean:
		rm *.o ; \
		rm irregular; 


The example assumes that you have the latest (per july 2022) version of the Basilisk library properly installed, and the `CFDwavemaker.h` and `libCFDwavemaker_openmp.a` file copied to the work directory. 

 

Linking up  with OpenFOAM
-------------------------

This is a bit more tricky since OpenFOAM comes in various versions and distributions. I have chosen the .com version of openFOAM, and developed some extensions to the built-in waveModel of OpenFOAM. The scripts will be available soon, in a separate repo.  








