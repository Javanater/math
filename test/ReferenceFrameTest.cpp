//
// Created by Madison on 8/26/2016.
//

#include <math/geometry/ReferenceFrame.hpp>
#include "gtest/gtest.h"

using namespace std;
using namespace flabs;
using namespace Eigen;

TEST(ReferenceFrameTest, 2d_unit_circle)
{
	ReferenceFrame<2> world;
	ReferenceFrame<2> point(1, 0, 0, &world);
	for (double       i = -M_PI; i < M_PI; i += .01)
	{
		world.setYawOffset(i);
		double x, y, t;
		point.getXYYaw(x, y, t);
		ASSERT_NEAR(cos(i), x, numeric_limits<double>::epsilon());
		ASSERT_NEAR(sin(i), y, numeric_limits<double>::epsilon());
		ASSERT_NEAR(i, t, numeric_limits<double>::epsilon());
	}
}
