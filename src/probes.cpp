#include "probes.h"


void Probes::init_probes() {

	if (dirExists( ts_path.c_str()) == 0) {
		std::cout << "WARNING: Specified directory for storage of time-series files does not exist. Directory will be created at the following path:  " << ts_path << std::endl;
		createDirectory(ts_path);
	}

	for (int i = 0; i < ts_nprobes; i++) {
		char buffer[256]; sprintf(buffer, "%05d", i);
		std::string str(buffer);
		std::string fpath = (ts_path + ts_filename + buffer + ".dat");
		//std::cout << fpath << std::endl;
		FILE* fp = fopen(fpath.c_str(), "w");
		// create file, append coordinates to top line
		fprintf(fp, "# CFDwavemaker probe %d, x= %.3f , y= %.3f , z= %.3f \n", i, coords[i].x, coords[i].y, coords[i].z);
		fprintf(fp, "# %12s  %14s  %14s  %14s  %14s\n", "Time", "wave_elev", "u", "v", "w");
		fclose(fp);
	}
}

bool Probes::checkTime(double tpt) {
	if (tpt >= update_time) {
		update_time += dt;
		return true;
	}
	return false;

}

void Probes::write(double tpt, Irregular& irregular, Ramp& ramp) {
	// loop over probes and write
	for (int i = 0; i < ts_nprobes; i++) {
		double xpt = coords[i].x;
		double ypt = coords[i].y;
		double zpt = coords[i].z;
		double welev, ux, uy, uz;

		welev = ramp.ramp(tpt, xpt, ypt) * irregular.eta(tpt, xpt, ypt);
		ux = ramp.ramp(tpt, xpt, ypt) * irregular.u(tpt, xpt, ypt, zpt);
		uy = ramp.ramp(tpt, xpt, ypt) * irregular.v(tpt, xpt, ypt, zpt);
		uz = ramp.ramp(tpt, xpt, ypt) * irregular.w(tpt, xpt, ypt, zpt);
						

		// open file and write line
		char buffer[256]; sprintf(buffer, "%05d", i);
		std::string str(buffer);
		std::string fpath = (ts_path + ts_filename + buffer + ".dat");
		//std::cout << fpath << std::endl;
		FILE* fp = fopen(fpath.c_str(), "a");
		fprintf(fp, "%14.4f  %14.4f  %14.4f  %14.4f  %14.4f\n", tpt, welev, ux, uy, uz);
		// create file, do nothing.
		fclose(fp);
	}
}

void Probes::write(double tpt, lsGrid& sgrid, Irregular& irregular, Ramp& ramp) {
	if (!sgrid.CheckTime(tpt)) {
#pragma omp single nowait
		sgrid.update(irregular, tpt);
	}

	// loop over probes and write
	for (int i = 0; i < ts_nprobes; i++) {
		double xpt = coords[i].x;
		double ypt = coords[i].y;
		double zpt = coords[i].z;
		double welev, ux, uy, uz;
		welev = ramp.ramp(tpt, xpt, ypt) * sgrid.eta(tpt, xpt, ypt);
		ux = ramp.ramp(tpt, xpt, ypt) * sgrid.u(tpt, xpt, ypt, zpt);
		uy = ramp.ramp(tpt, xpt, ypt) * sgrid.v(tpt, xpt, ypt, zpt);
		uz = ramp.ramp(tpt, xpt, ypt) * sgrid.w(tpt, xpt, ypt, zpt);
		std::cout << welev << std::endl;

		// open file and write line
		char buffer[256]; sprintf(buffer, "%05d", i);
		std::string str(buffer);
		std::string fpath = (ts_path + ts_filename + buffer + ".dat");
		//std::cout << fpath << std::endl;
		FILE* fp = fopen(fpath.c_str(), "a");
		fprintf(fp, "%14.4f  %14.4f  %14.4f  %14.4f  %14.4f\n", tpt, welev, ux, uy, uz);
		// create file, do nothing.
		fclose(fp);
	}
}

void Probes::write(double tpt, lsGridSpline& sgrids, Irregular& irregular, Ramp& ramp) {
	if (!sgrids.CheckTime(tpt)) {
#pragma omp single nowait
		sgrids.update(irregular, tpt);
	}

	// loop over probes and write
	for (int i = 0; i < ts_nprobes; i++) {
		std::vector<double> v;
		double xpt = coords[i].x;
		double ypt = coords[i].y;
		double zpt = coords[i].z;
		v = sgrids.get_kinematics_at_point(tpt, xpt, ypt, zpt, irregular.depth);
	
		
		double welev, ux, uy, uz;
		welev = ramp.ramp(tpt, xpt, ypt) * v[0];
		ux = ramp.ramp(tpt, xpt, ypt) * v[1];
		uy = ramp.ramp(tpt, xpt, ypt) * v[2];
		uz = ramp.ramp(tpt, xpt, ypt) * v[3];
		//std::cout << welev << std::endl;

		// open file and write line
		char buffer[256]; sprintf(buffer, "%05d", i);
		std::string str(buffer);
		std::string fpath = (ts_path + ts_filename + buffer + ".dat");
		//std::cout << fpath << std::endl;
		FILE* fp = fopen(fpath.c_str(), "a");
		fprintf(fp, "%14.4f  %14.4f  %14.4f  %14.4f  %14.4f\n", tpt, welev, ux, uy, uz);
		// create file, do nothing.
		fclose(fp);
	}
}

void Probes::write(double tpt, Wavemaker& wavemaker, Ramp& ramp){
	// loop over probes and write
	for (int i = 0; i < ts_nprobes; i++) {
		double xpt = coords[i].x;
		double ypt = coords[i].y;
		double zpt = coords[i].z;
		double welev, ux, uy, uz;
		welev = ramp.ramp(tpt, xpt, ypt) * wavemaker.wave_elev_piston(tpt);
		ux = ramp.ramp(tpt, xpt, ypt) * wavemaker.u_piston(tpt);
		uy = 0.;
		uz = 0.;

		// open file and write line
		char buffer[256]; sprintf(buffer, "%05d", i);
		std::string str(buffer);
		std::string fpath = (ts_path + ts_filename + buffer + ".dat");
		//std::cout << fpath << std::endl;
		FILE* fp = fopen(fpath.c_str(), "a");
		fprintf(fp, "%14.4f  %14.4f  %14.4f  %14.4f  %14.4f\n", tpt, welev, ux, uy, uz);
		// create file, do nothing.
		fclose(fp);
	}
}


void Probes::write(double tpt, Stokes5& stokes5, Ramp& ramp) {
	// loop over probes and write
	for (int i = 0; i < ts_nprobes; i++) {
		double xpt = coords[i].x;
		double ypt = coords[i].y;
		double zpt = coords[i].z;
		double welev, ux, uy, uz;
		welev = ramp.ramp(tpt, xpt, ypt) * stokes5.eta(tpt, xpt, ypt);
		ux = ramp.ramp(tpt, xpt, ypt) * stokes5.u(tpt, xpt, ypt, zpt);
		uy = ramp.ramp(tpt, xpt, ypt) * stokes5.v(tpt, xpt, ypt, zpt);
		uz = ramp.ramp(tpt, xpt, ypt) * stokes5.w(tpt, xpt, ypt, zpt);

		// open file and write line
		char buffer[256]; sprintf(buffer, "%05d", i);
		std::string str(buffer);
		std::string fpath = (ts_path + ts_filename + buffer + ".dat");
		//std::cout << fpath << std::endl;
		FILE* fp = fopen(fpath.c_str(), "a");
		fprintf(fp, "%14.4f  %14.4f  %14.4f  %14.4f  %14.4f\n", tpt, welev, ux, uy, uz);
		// create file, do nothing.
		fclose(fp);
	}
}


#if defined(SWD_enable)

void Probes::write(double tpt, SpectralWaveData* swd, Ramp& ramp){
	// Tell the swd object current application time...
	try {
		swd->UpdateTime(tpt);
	}
	catch (SwdInputValueException& e) {  //Could be t > tmax from file.
		std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
		// If we will try again with a new value of t
		// we first need to call: swd.ExceptionClear()
		exit(EXIT_FAILURE);  // In this case we just abort.
	}
	// loop over probes and write
	for (int i = 0; i < ts_nprobes; i++) {
		double xpt = coords[i].x;
		double ypt = coords[i].y;
		double zpt = coords[i].z;
		double welev, ux, uy, uz;
		//std::cout << "time: " << tpt << std::endl;
		welev = ramp.ramp(tpt, xpt, ypt) * swd->Elev(xpt, ypt);
		vector_swd U = swd->GradPhi(xpt, ypt, zpt);
		ux = ramp.ramp(tpt, xpt, ypt) * U.x;
		uy = ramp.ramp(tpt, xpt, ypt) * U.y;
		uz = ramp.ramp(tpt, xpt, ypt) * U.z;

		// open file and write line
		char buffer[256]; sprintf(buffer, "%05d", i);
		std::string str(buffer);
		std::string fpath = (ts_path + ts_filename + buffer + ".dat");
		//std::cout << fpath << std::endl;
		FILE* fp = fopen(fpath.c_str(), "a");
		fprintf(fp, "%14.4f  %14.4f  %14.4f  %14.4f  %14.4f\n", tpt, welev, ux, uy, uz);
		// create file, do nothing.
		fclose(fp);
	}
}

void Probes::write(double tpt, lsGrid& sgrid, SpectralWaveData* swd, Ramp& ramp) {
	if (!sgrid.CheckTime(tpt)) {
#pragma omp single nowait
		sgrid.update(swd, tpt);
	}
	// loop over probes and write
	for (int i = 0; i < ts_nprobes; i++) {
		double xpt = coords[i].x;
		double ypt = coords[i].y;
		double zpt = coords[i].z;
		double welev, ux, uy, uz;
		welev = ramp.ramp(tpt, xpt, ypt) * sgrid.eta(tpt, xpt, ypt);
		ux = ramp.ramp(tpt, xpt, ypt) * sgrid.u(tpt, xpt, ypt, zpt);
		uy = ramp.ramp(tpt, xpt, ypt) * sgrid.v(tpt, xpt, ypt, zpt);
		uz = ramp.ramp(tpt, xpt, ypt) * sgrid.w(tpt, xpt, ypt, zpt);

		// open file and write line
		char buffer[256]; sprintf(buffer, "%05d", i);
		std::string str(buffer);
		std::string fpath = (ts_path + ts_filename + buffer + ".dat");
		//std::cout << fpath << std::endl;
		FILE* fp = fopen(fpath.c_str(), "a");
		fprintf(fp, "%14.4f  %14.4f  %14.4f  %14.4f  %14.4f\n", tpt, welev, ux, uy, uz);
		// create file, do nothing.
		fclose(fp);
	}
}
#endif
