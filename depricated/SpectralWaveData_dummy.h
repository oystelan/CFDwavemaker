#ifndef SWD_CPP_H_INCLUDED
#define SWD_CPP_H_INCLUDED

/*
This is a dummy class to include instead of SWD in case of compilation on windows. 

Oystein Lande

*/

// Floating point model for interfacing with spectral_wave_data
#ifdef SWD_API_FLOAT
typedef float real_swd;
#else
typedef double real_swd;
#endif

typedef struct {
    real_swd x;
    real_swd y;
    real_swd z;
} vector_swd;

// Specific Exception Classes...

class SwdException : public std::runtime_error
{
public:
    SwdException(const char* msg) : std::runtime_error(msg) { }
};

class SwdFileCantOpenException : public SwdException
{
public:
    SwdFileCantOpenException(const char* msg) : SwdException(msg) { }
};

class SwdFileBinaryException : public SwdException
{
public:
    SwdFileBinaryException(const char* msg) : SwdException(msg) { }
};

class SwdFileDataException : public SwdException
{
public:
    SwdFileDataException(const char* msg) : SwdException(msg) { }
};

class SwdInputValueException : public SwdException
{
public:
    SwdInputValueException(const char* msg) : SwdException(msg) { }
};

class SwdAllocateException : public SwdException
{
public:
    SwdAllocateException(const char* msg) : SwdException(msg) { }
};

class SpectralWaveData
{
    void* obj; // Wrapper to the C-object of the ocean wave model

public:
    SpectralWaveData(std::string file_swd, real_swd x0, real_swd y0,
        real_swd t0, real_swd beta, real_swd rho = 1025.0,
        int nsumx = -1, int nsumy = -1, int impl = 0, int ipol = 0,
        int norder = 0, bool dc_bias = false) {
    }
    
    /*
    file_swd:          Name of actual swd file

    x0, y0, t0, beta:  Relation between SWD and application coordinates.
                       beta in degree.

    rho:               Density of water(applied for pressure calculations)
    nsumx, nsumy       Number of spectral components to apply (<0: apply all)

    impl               Index to determine actual derived class
                       0 = Default
                      <0 = In-house and experimental implementations
                      >0 = Validated implementations available open software

    ipol               Index to request actual temporal interpolation scheme
                       0 = Default (C^2 continous scheme)
                       1 = C^1 continous
                       2 = C^3 continous

    norder             Expansion order to apply in kinematics for z>0
                       0 = Apply expansion order specified in swd file (default)
                       <0 = Apply exp(kj z)
                       >0 = Apply expansion order = norder

    dc_bias            Control application of zero-frequency bias present in SWD file
                       false = Suppress contribution from zero frequency amplitudes (default)
                       true  = Apply zero frequency amplitudes from SWD file.

    C++ exceptions the constructor may throw:  (should be catched in application)

    SwdException:               Base class for the following exceptions
    SwdFileCantOpenException:   Typical if SWD file is not an existing file
    SwdFileBinaryException:     SWD file does not apply float/little-endian
    SwdFileDataException:       Error during reading and checking data
    SwdInputValueException:     Input arguments for class methods are not sound
    SwdAllocateException:       Not able to allocate internal SWD storage

    /**/

    ~SpectralWaveData();

    // =================================================================
    // Methods for evaluating kinematics...
    // =================================================================

    // Apply current application time
    // Possible exceptions thrown: SwdFileDataException, SwdInputValueException
    void UpdateTime(real_swd time) {}

    // Calculate velocity potential
    // Possible exceptions thrown: None
    real_swd Phi(real_swd x, real_swd y, real_swd z) {
        return 0.;
    }

    // Calculate stream function
    // Possible exceptions thrown: None
    real_swd Stream(real_swd x, real_swd y, real_swd z) {
        return 0.;
    }

    // Calculate time derivative of wave potential (earth fixed observer)
    // Possible exceptions thrown: None
    real_swd DdtPhi(real_swd x, real_swd y, real_swd z) {
        return 0.;
    }

    // Calculate wave elevation
    real_swd Elev(real_swd x, real_swd y)
    {
        return 0.;
    }

    // Calculate particle velocity
    // Possible exceptions thrown: None
    vector_swd GradPhi(real_swd x, real_swd y, real_swd z) {
        vector_swd dd;
        dd.x = dd.y = dd.z = 0.;
        return dd;
    }


    // Complete Bernoulli pressure
    // Possible exceptions thrown: None
    real_swd Pressure(real_swd x, real_swd y, real_swd z) {
        return 0.;
    }


    // ===================================================================
//  Provide parameters from the swd-file:

// Extract the character parameter 'name' from object
    std::string GetChr(std::string const& name)
    {
        return name.c_str();
    }

    // Extract the int parameter 'name' from object
    int GetInt(std::string const& name)
    {
        return 0;
    }

    // Extract the real parameter 'name' from object
    real_swd GetReal(std::string const& name)
    {
        return 0.;
    }
    // ===================================================================

    // Clear error flag in C/Fortran implementation.
    // Only applied in case of advanced exception handling. (recovery)
    void SpectralWaveData::ExceptionClear(){}

};
#endif //SWD_CPP_H_INCLUDED
