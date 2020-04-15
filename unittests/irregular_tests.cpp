#include "pch.h"

#include "irregular.h"
struct IrregularTest : testing::Test {
	Irregular* irregular;
	
	IrregularTest() {
		irregular = new Irregular;
		// Parameters
		irregular->nfreq = 3;
		// Variables
		irregular->nfreq = 3;
		irregular->ndir = 1;
		irregular->extrapolation_met = 0;
		irregular->bandwidth = 1000;
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

	double t = 0.;
	double x = 0.;
	double y = 0.;
	//std::cout << "Wave elevation: " << irregular->eta(0.0, 0.0, 0.0) << std::endl;
	EXPECT_DOUBLE_EQ(6., irregular->eta(0.0, 0.0, 0.0));
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
