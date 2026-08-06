# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 14:58:16 2016

@author: oland
"""

from ctypes import*
import matplotlib.pylab as plt
import numpy as np
import os
import shutil

import time

def TicTocGenerator():
    # Generator that returns time differences
    ti = 0           # initial time
    tf = time.time() # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti # returns the time difference

TicToc = TicTocGenerator() # create an instance of the TicTocGen generator

# This will be the main function through which we define both tic() and toc()
def toc(tempBool=True):
    # Prints the time difference yielded by generator instance TicToc
    tempTimeInterval = next(TicToc)
    if tempBool:
        print( "Elapsed time: %f seconds.\n" %tempTimeInterval )

def tic():
    # Records a time in TicToc, marks the beginning of a time interval
    toc(False)


# give location of dll
#mydll = cdll.comflow
temp = os.path.abspath(__file__)
temp = os.path.realpath(temp)
temp = os.path.dirname(temp)

dlclose_func = cdll.LoadLibrary('').dlclose
dlclose_func.argtypes = [c_void_p]

def dlopen(fpath):
    mydll = cdll.LoadLibrary(fpath)
    return mydll

def is_loaded(lib):
    libp = os.path.abspath(lib)
    ret = os.system("lsof -p %d | grep %s > /dev/null" % (os.getpid(), libp))
    return (ret == 0)

def reload_lib(lib):
    handle = lib._handle
    name = lib._name
    del lib
    while is_loaded(name):
        dlclose_func(handle)   
        #libdl = CDLL("libdl.so")
        #libdl.dlclose(handle)
    return cdll.LoadLibrary(name)

# def dlclose(mydll):
#     handle = mydll._handle
#     del mydll
#     dlclose_func(handle)

def dlclose(lib):
    handle = lib._handle
    name = lib._name
    del lib
    while is_loaded(name):
        dlclose_func(handle)   
        #libdl = CDLL("libdl.so")
        #libdl.dlclose(handle)       

# and define the functions which CFDwavemaker.h tells us are in the library

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

def wave_Kinematics(mydll,t,x,y,z):
    aa = mydll.wave_Kinematics
    #aa.restype = c_double
    aa.restype = POINTER(c_double*4)
    aa.argtypes = [c_double,c_double,c_double,c_double]
    res = aa(c_double(x),c_double(y),c_double(z),c_double(t)).contents
    #print(res[0], res[1], res[2], res[3])
    return res[0], res[1], res[2], res[3]

def wave_update(mydll, t):
    aa = mydll.wave_force_update
    aa.restype = c_int
    aa.argtype = c_double
    return(aa(c_double(t)))

def init_dll(mydll):
    aa = mydll.Init
    aa.restype = c_int
    aa.argtypes = [POINTER(c_double),POINTER(c_double)]
    return aa(c_double(0),c_double(0))

def clean_up(mydll):
    aa = mydll.Cleanup
    aa.restype = c_int
    return aa()


# All done. lets test:

os.environ["OMP_NUM_THREADS"] = "12" # export OMP_NUM_THREADS=4

tic()
x = 0.
xx = np.arange(0,100,0.1)
y = 0.0
z = -10

# Step 1 - Lest compute the exact second order solution
shutil.copyfile("waveinput_orig.dat", 'waveinput.dat')
mydll1 = dlopen("../../builds/linux64/libCFDwavemaker_swd.so")  
init_dll(mydll1)

wave_elev1 = []
u1 = []
v1 = []
w1 = []

for t in xx:
    #print(wave_SurfElev(mydll,t,-200,0))
    wave_elev1.append(wave_SurfElev(mydll1,t,x,y))
    u1.append(wave_VeloX(mydll1,t,x,y,z))
    v1.append(wave_VeloY(mydll1,t,x,y,z))
    w1.append(wave_VeloZ(mydll1,t,x,y,z))
clean_up(mydll1)
dlclose(mydll1)
print("done step 1")
toc()
tic()
# Step 2 - now use boost version of irreg
mydll2 = dlopen("../../builds/linux64/libCFDwavemaker_swd.so")  
shutil.copyfile("waveinput_boost.dat", 'waveinput.dat')
init_dll(mydll2)

wave_elev2 = []
u2 = []
v2 = []
w2 = []

for t in xx:
    wave_update(mydll2, t)
    #print(wave_SurfElev(mydll,t,-200,0))
    wave_elev2.append(wave_SurfElev(mydll2,t,x,y))
    u2.append(wave_VeloX(mydll2,t,x,y,z))
    v2.append(wave_VeloY(mydll2,t,x,y,z))
    w2.append(wave_VeloZ(mydll2,t,x,y,z))
clean_up(mydll2)
dlclose(mydll2)
toc()

 




plt.plot(xx, wave_elev1, label="orig")
plt.plot(xx, wave_elev2,"--", label="boost")
#plt.plot(xx, wave_elev3,"--", label="spline")
plt.legend()
plt.grid(True)
plt.savefig("./result_eta.png")
plt.show()

#exit()
plt.clf()
plt.plot(xx,u1, label='u orig')
plt.plot(xx,v1, label='v orig')
plt.plot(xx,w1, label='w orig')
plt.plot(xx,u2,'--', label='u boost')
plt.plot(xx,v2,'--', label='v boost')
plt.plot(xx,w2,'--', label='w boost')
#plt.plot(xx,u3,':', label='u spline')
#plt.plot(xx,v3,':', label='v spline')
#plt.plot(xx,w3,':', label='w spline')
plt.legend()
plt.grid(True)
plt.savefig("./result_u.png")
plt.show()
# all done. remember to clean up after us.






