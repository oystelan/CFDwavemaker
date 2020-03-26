Installation
============

This pages describes how to install and get up and running with the CFDwavemaker C++ library, with different CFD solvers.

Linking up with ComFLOW
-----------------------

.. code:: C++

   fprintf(stderr,"A literal block directive explicitly marked as python code\n");

Linking up with Basilisk (or C)
-------------------------------

Linking up  with OpenFOAM
-------------------------

Linking up with Python
----------------------

Occationally it may be useful generate kinematics data for other reasons than pushing it into a CFD solver. In this case python 

.. code:: python

	from ctypes import*
	import matplotlib.pylab as plt
	import numpy as np

	# give location of dll (this is used in windows)
	# mydll = cdll.CFDwavemaker

	# And this is for Linux/unix
	mydll = cdll.LoadLibrary("./CFDwavemaker.so") 

	def velocityx(mydll,t,x,y,z):
	    aa = mydll.VelocityX
	    aa.restype = c_double
	    aa.argtypes = [c_int,c_int,c_int,c_double,c_double,c_double,c_double]
	    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t))

	def velocityy(mydll,t,x,y,z):
	    aa = mydll.VelocityY
	    aa.restype = c_double
	    aa.argtypes = [c_int, c_int, c_int, c_double, c_double, c_double, c_double]
	    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t))

	def velocityz(mydll,t,x,y,z):
	    aa = mydll.VelocityZ
	    aa.restype = c_double
	    aa.argtypes = [c_int,c_int,c_int,c_double,c_double,c_double,c_double]
	    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t))

	def waveelev(mydll,t,x,y):
	    aa = mydll.SurfaceElevation
	    aa.restype = c_double
	    aa.argtypes = [c_int,c_int,c_double,c_double,c_double]
	    return aa(c_int(0),c_int(0),c_double(x),c_double(y),c_double(t))

	def volfrac(mydll,x,y,z,t,delta):
	    aa = mydll.VolumeFraction
	    aa.restype = c_double
	    aa.argtypes = [c_double,c_double,c_double,c_double,c_double]
	    return aa(c_double(x),c_double(y),c_double(z),c_double(t),c_double(delta))

	def init_dll(mydll):
	    aa = mydll.Init
	    aa.restype = c_int
	    aa.argtypes = [POINTER(c_double),POINTER(c_double)]
	    return aa(c_double(0),c_double(0));

	def clean_up(mydll):
	    aa = mydll.Cleanup
	    aa.restype = c_int
	    return aa();


	# We start of by calling the init function. This will read the waveinput.dat file and perform necessary initialization
	print(init_dll(mydll))

	# We create a time vector from 0 to 100 sec, with 0.5sec spacing
	time = np.arange(0, 100, 0.5)

	# lets extract the surface elevation at x=y=0
	wave_elev = np.zeros(len(time))
	for wave,t in zip(wave_elev,time):
	    wave = waveelev(mydll,t,0.0,0.0) # x and y position set to 0.

	# or perhaps the horizontal particle velocity at position x = 0, y=0, z = -10.

	ux = np.zeros(len(time))
	for u,t in zip(ux,time):
	    u = velocityx(mydll,t,0.0,0.0,-11.0) # x and y position set to 0.

	# You get the idea.

	# Plotting some kinematics data
	plt.plot(time,wave)
	plt.plot(time,ux)
	plt.show()

	# now that we are done, we close down the link to the library before exiting
	clean_up(mydll)







