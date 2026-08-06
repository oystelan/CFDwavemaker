## -*- coding: utf-8 -*-
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
    res = aa(c_double(x),c_double(y),c_double(z),c_double(t),c_double(t)).contents
    #print(res[0], res[1], res[2], res[3])
    return res[0], res[1], res[2], res[3], res[4]

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
os.environ["OMP_NUM_THREADS"] = "2" # export OMP_NUM_THREADS=4
x = np.arange(-32,32.,1)
z = 0.
t = 0.
y = 0.


tic()
# Step 1 - Lest compute the exact second order solution
#shutil.copyfile("waveinput_stokes.dat", 'waveinput.dat')
mydll1 = dlopen("../../builds/linux64/libCFDwavemaker_swd.so")  
init_dll(mydll1)

wave_elev1 = []
u1 = []
v1 = []
w1 = []

for xx in x:
    
    wave_elev1.append(wave_SurfElev(mydll1,t,xx,y))
    print(xx,wave_elev1[-1])
    #u1.append(wave_VeloX(mydll1,t,xx,y,z))
    #v1.append(wave_VeloY(mydll1,t,xx,y,z))
    #w1.append(wave_VeloZ(mydll1,t,xx,y,z))
clean_up(mydll1)
dlclose(mydll1)

toc()
# compute slope
#print("slope angle x: ", np.arctan(wave_elev1[-1]/x[-1])*180/np.pi)


# plt.figure(figsize=(20,6))
plt.plot(x, wave_elev1, label="cfdwavemaker")
#plt.plot(x, u1, x, v1, x, w1)
# plt.plot(time, wave_elev2,"--", label="raschii")
# plt.legend()
plt.grid(True)
# plt.savefig("./result_eta.png")
plt.show()






