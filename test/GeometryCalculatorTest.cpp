//
// Created by Madison on 8/26/2016.
//

#include <math/geometry/GeometryCalculator.hpp>
#include "gtest/gtest.h"

#define TEST_COUNT 100000

using namespace std;
using namespace flabs;
using namespace Eigen;

TEST(GeometryCalculatorTest, 2d_orthogonal)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d vector           = Vector2d::Random();
		Vector2d orthogonalVector = orthogonal(vector);
		ASSERT_EQ(0, vector.dot(orthogonalVector));
	}
}

TEST(GeometryCalculatorTest, 3d_orthogonal)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector3d vector           = Vector3d::Random();
		Vector3d orthogonalVector = orthogonal(vector);
		ASSERT_LE(vector.dot(orthogonalVector),
			numeric_limits<double>::epsilon() * 2);
	}
}

TEST(GeometryCalculatorTest, 4d_orthogonal)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector4d vector           = Vector4d::Random();
		Vector4d orthogonalVector = orthogonal(vector);
		ASSERT_LE(vector.dot(orthogonalVector),
			numeric_limits<double>::epsilon() * 2);
	}
}

TEST(GeometryCalculatorTest, 2d_Parallel)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d direction = Vector2d::Random();
		Vector2d p1        = Vector2d::Random();
		Vector2d p2        = Vector2d::Random();
		ASSERT_EQ(NONE, intersects(p1, direction, p2, direction));
	}
}

TEST(GeometryCalculatorTest, 2d_Parallel_Opposite)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d direction         = Vector2d::Random();
		Vector2d oppositeDirection = -direction;
		Vector2d p1                = Vector2d::Random();
		Vector2d p2                = Vector2d::Random();
		ASSERT_EQ(NONE, intersects(p1, direction, p2, oppositeDirection));
	}
}

TEST(GeometryCalculatorTest, 2d_Coincident)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d direction = Vector2d::Random();
		Vector2d position  = Vector2d::Random();
		ASSERT_EQ(COINCIDENT,
			intersects(position, direction, position, direction));
	}
}

TEST(GeometryCalculatorTest, 2d_Coincident_Opposite)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d direction         = Vector2d::Random();
		Vector2d oppositeDirection = -direction;
		Vector2d position          = Vector2d::Random();
		ASSERT_EQ(COINCIDENT,
			intersects(position, direction, position, oppositeDirection));
	}
}

TEST(GeometryCalculatorTest, 2d_Intersection)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        = intersectionPoint - position1;
		Vector2d position2         = Vector2d::Random();
		Vector2d direction2        = intersectionPoint - position2;
		ASSERT_EQ(INTERSECT,
			intersects(position1, direction1, position2, direction2));
	}
}

TEST(GeometryCalculatorTest, 2d_Intersection_One_Away)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        = position1 - intersectionPoint;
		Vector2d position2         = Vector2d::Random();
		Vector2d direction2        = intersectionPoint - position2;
		ASSERT_EQ(INTERSECT,
			intersects(position1, direction1, position2, direction2));
	}
}

TEST(GeometryCalculatorTest, 2d_Intersection_Two_Away)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        = intersectionPoint - position1;
		Vector2d position2         = Vector2d::Random();
		Vector2d direction2        = position2 - intersectionPoint;
		ASSERT_EQ(INTERSECT,
			intersects(position1, direction1, position2, direction2));
	}
}

TEST(GeometryCalculatorTest, 2d_Intersection_Both_Away)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        = position1 - intersectionPoint;
		Vector2d position2         = Vector2d::Random();
		Vector2d direction2        = position2 - intersectionPoint;
		ASSERT_EQ(INTERSECT,
			intersects(position1, direction1, position2, direction2));
	}
}

TEST(GeometryCalculatorTest, 2d_Distance)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        =
					 (intersectionPoint - position1).normalized();
		Vector2d position2         = Vector2d::Random();
		Vector2d direction2        = intersectionPoint - position2;
		double   actualDistance    = (intersectionPoint - position1).norm();
		double   dist;
		ASSERT_EQ(INTERSECT,
			distance(position1, direction1, position2, direction2, dist));
		ASSERT_NEAR(actualDistance, dist, actualDistance * .000000001);
	}
}

TEST(GeometryCalculatorTest, 2d_Distance_One_Away)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        =
					 (position1 - intersectionPoint).normalized();
		Vector2d position2         = Vector2d::Random();
		Vector2d direction2        = intersectionPoint - position2;
		double   actualDistance    = -(intersectionPoint - position1).norm();
		double   dist;
		ASSERT_EQ(INTERSECT,
			distance(position1, direction1, position2, direction2, dist));
		ASSERT_NEAR(actualDistance, dist,
			std::abs(actualDistance * .000000001));
	}
}

TEST(GeometryCalculatorTest, 2d_Distance_Two_Away)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        =
					 (intersectionPoint - position1).normalized();
		Vector2d position2         = Vector2d::Random();
		Vector2d direction2        = position2 - intersectionPoint;
		double   actualDistance    = (intersectionPoint - position1).norm();
		double   dist;
		ASSERT_EQ(INTERSECT,
			distance(position1, direction1, position2, direction2, dist));
		ASSERT_NEAR(actualDistance, dist, actualDistance * .00000001);
	}
}

TEST(GeometryCalculatorTest, 2d_Distance_Both_Away)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        =
					 (position1 - intersectionPoint).normalized();
		Vector2d position2         = Vector2d::Random();
		Vector2d direction2        = position2 - intersectionPoint;
		double   actualDistance    = -(intersectionPoint - position1).norm();
		double   dist;
		ASSERT_EQ(INTERSECT,
			distance(position1, direction1, position2, direction2, dist));
		ASSERT_NEAR(actualDistance, dist,
			std::abs(actualDistance * .000000001));
	}
}

TEST(GeometryCalculatorTest, 2d_Intersection_Equals)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        = intersectionPoint - position1;
		Vector2d position2         = Vector2d::Random();
		Vector2d direction2        = intersectionPoint - position2;
		Vector2d result;
		ASSERT_EQ(INTERSECT,
			intersection(position1, direction1, position2, direction2, result));
		ASSERT_LE((double) (result - intersectionPoint).norm(), 1e-10);
	}
}

TEST(GeometryCalculatorTest, 2d_Intersection_Equals_One_Away)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        = position1 - intersectionPoint;
		Vector2d position2         = Vector2d::Random();
		Vector2d direction2        = intersectionPoint - position2;
		Vector2d result;
		ASSERT_EQ(INTERSECT,
			intersection(position1, direction1, position2, direction2, result));
		ASSERT_LE((double) (result - intersectionPoint).norm(), 1e-10);
	}
}

TEST(GeometryCalculatorTest, 2d_Intersection_Equals_Two_Away)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        = intersectionPoint - position1;
		Vector2d position2         = Vector2d::Random();
		Vector2d direction2        = position2 - intersectionPoint;
		Vector2d result;
		ASSERT_EQ(INTERSECT,
			intersection(position1, direction1, position2, direction2, result));
		ASSERT_LE((double) (result - intersectionPoint).norm(), 1e-10);
	}
}

TEST(GeometryCalculatorTest, 2d_Intersection_Equals_Both_Away)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        = position1 - intersectionPoint;
		Vector2d position2         = Vector2d::Random();
		Vector2d direction2        = position2 - intersectionPoint;
		Vector2d result;
		ASSERT_EQ(INTERSECT,
			intersection(position1, direction1, position2, direction2, result));
		ASSERT_LE((double) (result - intersectionPoint).norm(), 1e-10);
	}
}

TEST(GeometryCalculatorTest, 2d_Intersection_Equals_d1_d2)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        =
					 (intersectionPoint - position1).normalized();
		Vector2d position2         = Vector2d::Random();
		Vector2d direction2        =
					 (intersectionPoint - position2).normalized();
		double   ad1               = (intersectionPoint - position1).norm();
		double   ad2               = (intersectionPoint - position2).norm();
		double   d1;
		double   d2;
		Vector2d result;
		ASSERT_EQ(INTERSECT,
			intersectionDistance(position1, direction1, position2, direction2,
				result, d1, d2));
		ASSERT_LE((double) (result - intersectionPoint).norm(), 1e-10);
		ASSERT_NEAR(ad1, d1, std::abs(ad1 * .000000001));
		ASSERT_NEAR(ad2, d2, std::abs(ad2 * .000000001));
	}
}

TEST(GeometryCalculatorTest, 2d_Intersection_Equals_d1_d2_One_Away)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        =
					 (position1 - intersectionPoint).normalized();
		Vector2d position2         = Vector2d::Random();
		Vector2d direction2        =
					 (intersectionPoint - position2).normalized();
		double   ad1               = -(intersectionPoint - position1).norm();
		double   ad2               = (intersectionPoint - position2).norm();
		double   d1;
		double   d2;
		Vector2d result;
		ASSERT_EQ(INTERSECT,
			intersectionDistance(position1, direction1, position2, direction2,
				result, d1, d2));
		ASSERT_LE((double) (result - intersectionPoint).norm(), 1e-10);
		ASSERT_NEAR(ad1, d1, std::abs(ad1 * .000000001));
		ASSERT_NEAR(ad2, d2, std::abs(ad2 * .000000001));
	}
}

TEST(GeometryCalculatorTest, 2d_Intersection_Equals_d1_d2_Two_Away)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        =
					 (intersectionPoint - position1).normalized();
		Vector2d position2         = Vector2d::Random();
		Vector2d direction2        =
					 (position2 - intersectionPoint).normalized();
		double   ad1               = (intersectionPoint - position1).norm();
		double   ad2               = -(intersectionPoint - position2).norm();
		double   d1;
		double   d2;
		Vector2d result;
		ASSERT_EQ(INTERSECT,
			intersectionDistance(position1, direction1, position2, direction2,
				result, d1, d2));
		ASSERT_LE((double) (result - intersectionPoint).norm(), 1e-10);
		ASSERT_NEAR(ad1, d1, std::abs(ad1 * .000000001));
		ASSERT_NEAR(ad2, d2, std::abs(ad2 * .000000001));
	}
}

TEST(GeometryCalculatorTest, 2d_Intersection_Equals_d1_d2_Both_Away)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        =
					 (position1 - intersectionPoint).normalized();
		Vector2d position2         = Vector2d::Random();
		Vector2d direction2        =
					 (position2 - intersectionPoint).normalized();
		double   ad1               = -(intersectionPoint - position1).norm();
		double   ad2               = -(intersectionPoint - position2).norm();
		double   d1;
		double   d2;
		Vector2d result;
		ASSERT_EQ(INTERSECT,
			intersectionDistance(position1, direction1, position2, direction2,
				result, d1, d2));
		ASSERT_LE((double) (result - intersectionPoint).norm(), 1e-10);
		ASSERT_NEAR(ad1, d1, std::abs(ad1 * .000000001));
		ASSERT_NEAR(ad2, d2, std::abs(ad2 * .000000001));
	}
}

TEST(GeometryCalculatorTest, 2d_Ray_Seg_Hit)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        =
					 (intersectionPoint - position1).normalized();
		Vector2d position2         = Vector2d::Random();
		Vector2d position3         = (intersectionPoint - position2) * 1.01;
		double   ad1               = (intersectionPoint - position1).norm();
		double   d1;
		double   d2;
		Vector2d result;
		ASSERT_EQ(INTERSECT,
			intersectionDistance(position1, direction1, position2, position3,
				result, d1, d2));
		ASSERT_LE((double) (result - intersectionPoint).norm(), 1e-10);
		ASSERT_NEAR(ad1, d1, std::abs(ad1 * .000000001));
		ASSERT_GT(1, d2);
		ASSERT_LT(0, d2);
	}
}

TEST(GeometryCalculatorTest, 2d_Ray_Seg_Hit_Behind)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        =
					 (position1 - intersectionPoint).normalized();
		Vector2d position2         = Vector2d::Random();
		Vector2d position3         = (intersectionPoint - position2) * 1.01;
		double   ad1               = -(intersectionPoint - position1).norm();
		double   d1;
		double   d2;
		Vector2d result;
		ASSERT_EQ(INTERSECT,
			intersectionDistance(position1, direction1, position2, position3,
				result, d1, d2));
		ASSERT_LE((double) (result - intersectionPoint).norm(), 1e-10);
		ASSERT_NEAR(ad1, d1, std::abs(ad1 * .000000001));
		ASSERT_GT(1, d2);
		ASSERT_LT(0, d2);
	}
}

TEST(GeometryCalculatorTest, 2d_Ray_Seg_Miss)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        =
					 (intersectionPoint - position1).normalized();
		Vector2d position2         = Vector2d::Random();
		Vector2d position3         = (intersectionPoint - position2) * .99;
		double   ad1               = (intersectionPoint - position1).norm();
		double   d1;
		double   d2;
		Vector2d result;
		ASSERT_EQ(INTERSECT,
			intersectionDistance(position1, direction1, position2, position3,
				result, d1, d2));
		ASSERT_LE((double) (result - intersectionPoint).norm(), 1e-10);
		ASSERT_NEAR(ad1, d1, std::abs(ad1 * .000000001));
		ASSERT_LT(1, d2);
	}
}

TEST(GeometryCalculatorTest, 2d_Ray_Seg_Miss_Behind)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        =
					 (intersectionPoint - position1).normalized();
		Vector2d position2         = Vector2d::Random();
		Vector2d position3         = (intersectionPoint - position2) * 1.01;
		position2 += position3;
		double   ad1 = (intersectionPoint - position1).norm();
		double   d1;
		double   d2;
		Vector2d result;
		ASSERT_EQ(INTERSECT,
			intersectionDistance(position1, direction1, position2, position3,
				result, d1, d2));
		ASSERT_LE((double) (result - intersectionPoint).norm(), 1e-10);
		ASSERT_NEAR(ad1, d1, std::abs(ad1 * .000000001));
		ASSERT_GT(0, d2);
	}
}

TEST(GeometryCalculatorTest, 2d_Ray_Seg_Behind_Miss)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        =
					 (position1 - intersectionPoint).normalized();
		Vector2d position2         = Vector2d::Random();
		Vector2d position3         = (intersectionPoint - position2) * .99;
		double   ad1               = -(intersectionPoint - position1).norm();
		double   d1;
		double   d2;
		Vector2d result;
		ASSERT_EQ(INTERSECT,
			intersectionDistance(position1, direction1, position2, position3,
				result, d1, d2));
		ASSERT_LE((double) (result - intersectionPoint).norm(), 1e-10);
		ASSERT_NEAR(ad1, d1, std::abs(ad1 * .000000001));
		ASSERT_LT(1, d2);
	}
}

TEST(GeometryCalculatorTest, 2d_Ray_Seg_Hehind_Miss_Behind)
{
	for (int i = 0; i < TEST_COUNT; ++i)
	{
		Vector2d intersectionPoint = Vector2d::Random();
		Vector2d position1         = Vector2d::Random();
		Vector2d direction1        =
					 (position1 - intersectionPoint).normalized();
		Vector2d position2         = Vector2d::Random();
		Vector2d position3         = (intersectionPoint - position2) * 1.01;
		position2 += position3;
		double   ad1 = -(intersectionPoint - position1).norm();
		double   d1;
		double   d2;
		Vector2d result;
		ASSERT_EQ(INTERSECT,
			intersectionDistance(position1, direction1, position2, position3,
				result, d1, d2));
		ASSERT_LE((double) (result - intersectionPoint).norm(), 1e-10);
		ASSERT_NEAR(ad1, d1, std::abs(ad1 * .000000001));
		ASSERT_GT(0, d2);
	}
}
