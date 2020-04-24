#include "pch.h"

#include "irregular.h"
#include "sgrid.h"


struct LSgridTest : testing::Test {
	Irregular* irregular;
	sGrid* sgrid;
	LSgridTest(){
		irregular = new Irregular;
		sgrid = new sGrid;
		
		// Parameters
		irregular->nfreq = 3;
		// Variables
		irregular->nfreq = 3;
		irregular->ndir = 1;
		irregular->extrapolation_met = 0;
		irregular->normalize = 0;
		irregular->order = 2;

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
		irregular->A.push_back(3.1);
		irregular->A.push_back(2.1);
		irregular->A.push_back(1.1);
		irregular->phase.push_back(0.);
		irregular->phase.push_back(0.);
		irregular->phase.push_back(0.);
		irregular->theta.push_back(PI / 8.);
		irregular->theta.push_back(0.);
		irregular->theta.push_back(0.);

		
		sgrid->water_depth = irregular->depth;
		
	}
	~LSgridTest() {
		delete irregular;
		delete sgrid;
	}

};


/* Compare grid interpolation function LSgrid with direct calculation from the irregular class*/
TEST_F(LSgridTest, test1) {

	// Set time and location.
	double x = 0.;
	double y = 0.;
	double z = -10.;
	double t = 4.;
	sgrid->t0 = t;

	// Set grid refinement parameters
	sgrid->domain[0] = 0.;
	sgrid->domain[1] = 30.;
	sgrid->domain[2] = 0.;
	sgrid->domain[3] = 20.;
	sgrid->nx = 61;
	sgrid->ny = 11;
	sgrid->nl = 20;
	//sgrid->tan_a = 0.01;
	sgrid->tan_b = 2.0;
	sgrid->allocate();

	// update grids
	sgrid->initialize_surface_elevation(*irregular, sgrid->t0);
	sgrid->initialize_kinematics(*irregular);


	// Surface elevation at point
	//std::cout << "Wave elevation: " << irregular->eta(t, x, y) << ", sgrid: " << sgrid->eta(x, y) << std::endl;
	EXPECT_DOUBLE_EQ(irregular->eta(t, x, y), sgrid->eta(t, x, y));

	// Kinematics at the surface comparison
	//std::cout << "Kinematics at the free surface for position x=" << x << ", y=" << y << std::endl;
	//std::cout << "ux: " << irregular->u(t, x, y, irregular->eta(t, x, y)) << ", sgrid: " << sgrid->u(t, x, y, sgrid->eta(t, x, y)) << std::endl;
	//std::cout << "vx: " << irregular->v(t, x, y, irregular->eta(t, x, y)) << ", sgrid: " << sgrid->v(t, x, y, sgrid->eta(t, x, y)) << std::endl;
	//std::cout << "wx: " << irregular->w(t, x, y, irregular->eta(t, x, y)) << ", sgrid: " << sgrid->w(t, x, y, sgrid->eta(t, x, y)) << std::endl;
	EXPECT_NEAR(irregular->u(t, x, y, irregular->eta(t, x, y)), sgrid->u(t, x, y, sgrid->eta(t, x, y)), 1E-5);
	EXPECT_NEAR(irregular->v(t, x, y, irregular->eta(t, x, y)), sgrid->v(t, x, y, sgrid->eta(t, x, y)), 1E-5);
	EXPECT_NEAR(irregular->w(t, x, y, irregular->eta(t, x, y)), sgrid->w(t, x, y, sgrid->eta(t, x, y)), 1E-5);

	// Kinematics at seabed:
	//std::cout << "Kinematics at the seabed h=" << sgrid->water_depth << ", for position x=" << x << ", y=" << y << std::endl;
	//std::cout << "ux: " << irregular->u(t, x, y, -irregular->depth) << ", sgrid: " << sgrid->u(t, x, y, -sgrid->water_depth) << std::endl;
	//std::cout << "vx: " << irregular->v(t, x, y, -irregular->depth) << ", sgrid: " << sgrid->v(t, x, y, -sgrid->water_depth) << std::endl;
	//std::cout << "wx: " << irregular->w(t, x, y, -irregular->depth) << ", sgrid: " << sgrid->w(t, x, y, -sgrid->water_depth) << std::endl;
	EXPECT_NEAR(irregular->u(t, x, y, -irregular->depth), sgrid->u(t, x, y, -sgrid->water_depth), 1E-5);
	EXPECT_NEAR(irregular->v(t, x, y, -irregular->depth), sgrid->v(t, x, y, -sgrid->water_depth), 1E-5);
	EXPECT_NEAR(irregular->w(t, x, y, -irregular->depth), sgrid->w(t, x, y, -sgrid->water_depth), 1E-5);

	// Kinematics at position = z:
	//std::cout << "Kinematics at the position z=" << z << ", for position x=" << x << ", y=" << y << std::endl;
	std::cout << "ux: " << irregular->u(t, x, y, z) << ", sgrid: " << sgrid->u(t, x, y, z) << std::endl;
	std::cout << "vx: " << irregular->v(t, x, y, z) << ", sgrid: " << sgrid->v(t, x, y, z) << std::endl;
	std::cout << "wx: " << irregular->w(t, x, y, z) << ", sgrid: " << sgrid->w(t, x, y, z) << std::endl;
	EXPECT_NEAR(irregular->u(t, x, y, z), sgrid->u(t, x, y, z), 1E-1);
	EXPECT_NEAR(irregular->v(t, x, y, z), sgrid->v(t, x, y, z), 1E-1);
	EXPECT_NEAR(irregular->w(t, x, y, z), sgrid->w(t, x, y, z), 1E-1);

	EXPECT_TRUE(true);
}


TEST_F(LSgridTest, vtueksporter) {
	double t = 4.;
	sgrid->t0 = t;

	// Set grid refinement parameters
	sgrid->domain[0] = 0.;
	sgrid->domain[1] = 300.;
	sgrid->domain[2] = 0.;
	sgrid->domain[3] = 100.;
	sgrid->nx = 61;
	sgrid->ny = 31;
	sgrid->nl = 10;
	//sgrid->tan_a = 0.01;
	sgrid->tan_b = 1.5;
	sgrid->allocate();

	// update grids
	sgrid->initialize_surface_elevation(*irregular, sgrid->t0);
	sgrid->initialize_kinematics(*irregular);

	// Write vtu file
	//char name[100];
	//sprintf(name, "./unittest.vtu");
	FILE* fp = fopen("./unittest.vtu", "w");
	sgrid->export_vtu(fp);
	fclose(fp);
	EXPECT_TRUE(true);
}