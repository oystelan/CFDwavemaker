#include "pch.h"
#include "omp.h"

#include "irregular.h"
struct IrregularTest : testing::Test {
	Irregular* irregular;
	
	IrregularTest() {
		irregular = new Irregular;
		// Variables
		irregular->nfreq = 3;
		irregular->ndir = 1;
		irregular->extrapolation_met = 0;
		irregular->normalize = 0;
		irregular->order = 1;

		irregular->ampl = 1.;
		irregular->depth = 300.;
		irregular->mtheta = 10.;
		irregular->tofmax = 0.;
		irregular->fpoint[0] = 0.;
		irregular->fpoint[1] = 0.;
		irregular->swl = 0.;

		irregular->omega.push_back(2. * PI / 10.);
		irregular->omega.push_back(2. * PI / 5.);
		irregular->omega.push_back(2. * PI / 4.);
		irregular->k.push_back(std::pow(2. * PI / 10., 2) / 9.81);
		irregular->k.push_back(std::pow(2. * PI / 5., 2) / 9.81);
		irregular->k.push_back(std::pow(2. * PI / 4., 2) / 9.81);
		irregular->A.push_back(3.0);
		irregular->A.push_back(2.0);
		irregular->A.push_back(1.0);
		irregular->phase.push_back(0.);
		irregular->phase.push_back(0.);
		irregular->phase.push_back(0.);
		irregular->theta.push_back(PI / 8.);
		irregular->theta.push_back(0.);
		irregular->theta.push_back(0.);
	}
	~IrregularTest() {
		delete irregular;
	}

};
TEST_F(IrregularTest, test1) {
	// check that linear wave elevtion sums up correctly
	double t = 0.;
	double x = 0.;
	double y = 0.;

	//std::cout << "Wave elevation: " << irregular->eta(0.0, 0.0, 0.0) << std::endl;
	EXPECT_DOUBLE_EQ(6., irregular->eta(0.0, 0.0, 0.0));

	irregular->print();
	EXPECT_TRUE(true);
}


TEST_F(IrregularTest, test2) {
	// Check that normalize works

	irregular->normalize = 1;
	irregular->normalize_data();
	double t = 0.;
	double x = 0.;
	double y = 0.;
	//std::cout << "Wave elevation: " << irregular->eta(0.0, 0.0, 0.0) << std::endl;
	EXPECT_DOUBLE_EQ(1., irregular->eta(0.0, 0.0, 0.0));

	irregular->ampl = 3.;
	irregular->normalize_data();
	EXPECT_DOUBLE_EQ(3., irregular->eta(0.0, 0.0, 0.0));

	EXPECT_TRUE(true);
}


TEST_F(IrregularTest, interpolatetest) {
	// Original data
	std::vector<double> xData = { 1, 5, 10, 15, 20 };
	std::vector<double> yData = { 0.3, 0.5, 0.8, 0.1, 0.14 };

	// Set up some points for interpolation in xVals
	const int NPTS = 20;
	std::vector<double> xVals, yVals;
	for (int i = 1; i <= NPTS; i++) xVals.push_back((double)i);

	// Interpolate
	for (double x : xVals)
	{
		double y = irregular->interpolate(xData, yData, x, true);
		yVals.push_back(y);
	}

	// Output
	#define SP << std::fixed << std::setw( 15 ) << std::setprecision( 6 ) <<
	#define NL << '\n'
	std::cout << "Original data:\n";
	for (int i = 0; i < xData.size(); i++) std::cout SP xData[i] SP yData[i] NL;
	std::cout << "\nInterpolated data:\n";
	for (int i = 0; i < xVals.size(); i++) std::cout SP xVals[i] SP yVals[i] NL;

	

	EXPECT_TRUE(true);
}




