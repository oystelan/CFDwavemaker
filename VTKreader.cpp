#include "VTKreader.h"

// Function for listing vtk files in directory
vector<string>* VTKreader::listdir(const char* dirname, const char* suffix, int& IERR) {
	DIR* dp;
	dirent* d;
	vector<string>* vec = new vector<string>;

	dp = opendir(dirname);
	if (dp) {
		while ((d = readdir(dp)) != NULL) {

			// you want here. Something like:
			if (strstr(d->d_name, suffix)) {
				//cout << "found an .cxx file: " << d->d_name << "\n";
				vec->push_back(d->d_name);
			}
		}
		sort(vec->begin(), vec->end());
		IERR = 0;
	}
	else {
		IERR = -1;
	}
	return vec;
}

// Function loading data from the frist file into data set 2. This function is only called during initialization
void VTKreader::loadInit(string path, const char* fname) {
	// Load file 1


	path.append(fname);
	cout << path << endl;

	//reader1 = vtkXMLStructuredGridReader::New();
	reader1 = vtkSmartPointer<vtkXMLStructuredGridReader>::New();
	reader1->SetFileName(path.c_str());
	reader1->Update();
	dataset1 = reader1->GetOutput();
	dataset1->GetBounds(bounds);
	//cout << "loaded successfully!\n";

	cout << "The bounds: " << bounds[0] << " " << bounds[1] << " " << bounds[2] << " " << bounds[3] << " " << bounds[4] << " " << bounds[5] << " " << endl;

	if (bounds[2] == bounds[3]) {
		cout << "2D modus switched on." << endl;
		switch2d = 1;
	}

}

// Function moving data from dataset2 to dataset1, and loads new step into dataset2
void VTKreader::loadNext(string path, const char* fname) {

	vtkDataArray* data;
	// Load file 1
	path.append(fname);
	//cout << path << endl;

	reader0 = reader1;
	reader0->Update();

	dataset0 = reader0->GetOutput(); // move dataset1 to dataset0	
					 //dataset0 = dataset1;
	//reader1 = vtkXMLStructuredGridReader::New();
	reader1 = vtkSmartPointer<vtkXMLStructuredGridReader>::New();
	reader1->SetFileName(path.c_str());
	reader1->Update();

	dataset1 = reader1->GetOutput(); // load new dataset into dataset1
									 //cout << "loaded successfully!\n";

									 // load time variable from the datasets
									 // Extract Time from dataset of file 1
	data = dataset0->GetFieldData()->GetArray("TIME");
	doubledata = vtkDoubleArray::SafeDownCast(data);
	time0 = doubledata->GetValue(0);
	data = dataset1->GetFieldData()->GetArray("TIME");
	doubledata = vtkDoubleArray::SafeDownCast(data);
	time1 = doubledata->GetValue(0);
	cout << "Time interval: " << time0 << " to " << time1 << " sec" << endl;
	//data->Delete();
	//cstring vfraqstr(chF.c_str());
	//cstring velostr(chV.c_str());
	// Check that the later timestep file contains the field variables needed


	vtkIdType numberOfPointArrays = dataset1->GetPointData()->GetNumberOfArrays();
	vfraq_field_located = 0;
	velo_field_located = 0;
	for (vtkIdType i = 0; i < numberOfPointArrays; i++)
	{

		if (strcmp(dataset1->GetPointData()->GetArrayName(i), Uname.c_str()) == 0) {
			velo_field_located = 1;
		}
	}
}
