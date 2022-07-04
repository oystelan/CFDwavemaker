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

A simple example is shown below, where


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


Static linking
..............


Dynamic linking
...............

 

Linking up  with OpenFOAM
-------------------------









