//
// Created by Madison on 8/26/2016.
//

#include <math/Math.hpp>
#include "gtest/gtest.h"

#define TEST_COUNT 100000

using namespace std;
using namespace flabs;
using namespace boost;
using namespace boost::math::constants;
using flabs::hypot;

TEST(MathTest, angleDifference)
{
	ASSERT_EQ(0, angleDifference(0, 0));
	ASSERT_EQ(pi<double>(), angleDifference(pi<double>(), 0.));
	ASSERT_EQ(0, angleDifference(two_pi<double>(), 0.));
	ASSERT_EQ(0, angleDifference(2 * two_pi<double>(), 0.));
	ASSERT_NEAR(0.1, angleDifference(0., two_pi<double>() - .1), 1e-9);
	ASSERT_NEAR(0.1, angleDifference(0.2, 0.1), 1e-9);
	ASSERT_NEAR(0.1, angleDifference(0.2 + two_pi<double>(), 0.1), 1e-9);
	ASSERT_NEAR(0.1, angleDifference(0.2 - two_pi<double>(), 0.1), 1e-9);
	ASSERT_NEAR(0.1, angleDifference(0.2, 0.1 + two_pi<double>()), 1e-9);
	ASSERT_NEAR(0.1, angleDifference(0.2, 0.1 - two_pi<double>()), 1e-9);
	ASSERT_NEAR(-0.1, angleDifference(0.1, 0.2), 1e-9);
	ASSERT_NEAR(-0.1, angleDifference(0.1 - two_pi<double>(), 0.2), 1e-9);
	ASSERT_NEAR(-0.1, angleDifference(0.1 + two_pi<double>(), 0.2), 1e-9);
	ASSERT_NEAR(-0.1, angleDifference(0.1, 0.2 + two_pi<double>()), 1e-9);
	ASSERT_NEAR(-0.1, angleDifference(0.1, 0.2 - two_pi<double>()), 1e-9);
}

TEST(MathTest, intervalDifference)
{
	double val = intervalDifference(0., 0., -pi<double>(), pi<double>());
	ASSERT_EQ(0, val);
	ASSERT_EQ(-pi<double>(),
		intervalDifference(pi<double>(), 0., -pi<double>(), pi<double>()));
	ASSERT_EQ(0,
		intervalDifference(two_pi<double>(), 0., -pi<double>(), pi<double>()));
	ASSERT_EQ(0, intervalDifference(2 * two_pi<double>(), 0., -pi<double>(),
		pi<double>()));
	ASSERT_NEAR(0.1, intervalDifference(0.2, 0.1, -pi<double>(), pi<double>()),
		1e-9);
	ASSERT_NEAR(0.1,
		intervalDifference(0.2 + two_pi<double>(), 0.1, -pi<double>(),
			pi<double>()), 1e-9);
	ASSERT_NEAR(0.1,
		intervalDifference(0.2 - two_pi<double>(), 0.1, -pi<double>(),
			pi<double>()), 1e-9);
	ASSERT_NEAR(0.1,
		intervalDifference(0.2, 0.1 + two_pi<double>(), -pi<double>(),
			pi<double>()), 1e-9);
	ASSERT_NEAR(0.1,
		intervalDifference(0.2, 0.1 - two_pi<double>(), -pi<double>(),
			pi<double>()), 1e-9);
	ASSERT_NEAR(-0.1, intervalDifference(0.1, 0.2, -pi<double>(), pi<double>()),
		1e-9);
	ASSERT_NEAR(-0.1,
		intervalDifference(0.1 - two_pi<double>(), 0.2, -pi<double>(),
			pi<double>()), 1e-9);
	ASSERT_NEAR(-0.1,
		intervalDifference(0.1 + two_pi<double>(), 0.2, -pi<double>(),
			pi<double>()), 1e-9);
	ASSERT_NEAR(-0.1,
		intervalDifference(0.1, 0.2 + two_pi<double>(), -pi<double>(),
			pi<double>()), 1e-9);
	ASSERT_NEAR(-0.1,
		intervalDifference(0.1, 0.2 - two_pi<double>(), -pi<double>(),
			pi<double>()), 1e-9);
}

TEST(MathTest, hypot)
{
	for (double a = -1000; a <= 1000; a += .5)
		ASSERT_NEAR(sqrt(a * a), hypot(a), 1e-9);

	for (double a = -1000; a <= 1000; a += .5)
		for (double b = -1000; b <= 1000; b += .5)
			ASSERT_NEAR(sqrt(a * a + b * b), hypot(a, b), 1e-9);

	for (double a = -100; a <= 100; a += .5)
		for (double b = -100; b <= 100; b += .5)
			for (double c = -100; c <= 100; c += .5)
				ASSERT_NEAR(sqrt(a * a + b * b + c * c), hypot(a, b, c), 1e-9);

	for (double a = -10; a <= 10; a += .5)
		for (double b = -10; b <= 10; b += .5)
			for (double c = -10; c <= 10; c += .5)
				for (double d = -10; d <= 10; d += .5)
					ASSERT_NEAR(sqrt(a * a + b * b + c * c + d * d),
						hypot(a, b, c, d), 1e-9);
}

TEST(MathTest, zeroTo2Pi)
{
	for (double a = 0; a < two_pi<double>(); a += .1)
		for (int i = -10; i <= 10; ++i)
			ASSERT_NEAR(a, zeroTo2Pi(a + two_pi<double>() * i), 1e-9);
}
