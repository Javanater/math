/*
 * GeometryCalculator.h
 *
 *  Created on: Mar 12, 2016
 *      Author: Madison
 */

#ifndef GEOMETRYCALCULATOR_H_
#define GEOMETRYCALCULATOR_H_

#include <limits>
#include <Eigen/Eigen>
#include <ostream>
#include <iostream>

namespace flabs
{
	enum IntersectionType
	{
		INTERSECT, COINCIDENT, NONE
	};

	std::ostream&
	operator<<(std::ostream& output, const IntersectionType& type);

	template<class T>
	inline T orthogonal2d(T v)
	{
		std::swap(v(0, 0), v(1, 0));
		v(0, 0) = -v(0, 0);
		return v;
	}

	template<class T>
	inline T orthogonal3d(const T& v)
	{
		T normal[2];
		normal[0](0, 0) = -v(1, 0) - v(2, 0);
		normal[0](1, 0) = v(0, 0);
		normal[0](2, 0) = v(0, 0);
		normal[1](0, 0) = -v(0, 0) - v(2, 0);
		normal[1](1, 0) = v(1, 0);
		normal[1](2, 0) = v(1, 0);
		int select = v(1, 0) == 0 && v(0, 0) == -v(2, 0);
		return normal[select];
	}

	//TODO: Templatize, so all dimensions are optimized at compile time.
	template<class T, class Scaler = typename T::Scalar>
	inline T orthogonal(const T& v,
		Scaler tolerance = std::numeric_limits<Scaler>::epsilon() * 4)
	{
		if (T::RowsAtCompileTime == 2)
			return orthogonal2d(v);
		else if (T::RowsAtCompileTime == 3)
			return orthogonal3d(v);
		else
		{
			T normal;

			for (int skip = 0; skip < T::RowsAtCompileTime; ++skip)
			{
				Scaler   sum = 0;
				for (int i   = 0; i < T::RowsAtCompileTime; ++i)
				{
					if (i != skip)
					{
						normal(i, 0) = v(skip, 0);
						sum -= v(i, 0);
					}
				}
				normal(skip, 0) = sum;
				if (sum != 0 && v(skip, 0) != 0)
					break;
			}
			return normal;
		}
	}

	/**
	 * Calculates the intersection type between 2 infinitely long, bidirectional
	 * lines.
	 */
	template<class T, class Scaler = typename T::Scalar>
	IntersectionType
	intersects(const T& p1, const T& v1, const T& p2, const T& v2,
		Scaler tolerance = std::numeric_limits<Scaler>::epsilon() * 4)
	{
		T      n2          = orthogonal(v2);
		Scaler denominator = v1.dot(n2);

		if (std::abs(denominator) <= tolerance)
		{
			Scaler numerator = (p2 - p1).dot(n2);
			if (std::abs(numerator) <= tolerance)
				return COINCIDENT;
			else
				return NONE;
		}

		return INTERSECT;
	}

	/**
	 * Calculates the distance from p1 along the first line to the second line.
	 * This is useful for when you do not care about the intersection point,
	 * because it is never calculated. The distance is in units of ||v1||.
	 */
	template<class T, class Scaler = typename T::Scalar>
	IntersectionType
	distance(const T& p1, const T& v1, const T& p2, const T& v2, Scaler& result,
		Scaler tolerance = std::numeric_limits<Scaler>::epsilon() * 4)
	{
		T      n2          = orthogonal(v2);
		Scaler denominator = v1.dot(n2);
		Scaler numerator   = (p2 - p1).dot(n2);

		if (std::abs(denominator) <= tolerance)
		{
			if (std::abs(numerator) <= tolerance)
			{
				result = 0;
				return COINCIDENT;
			}
			else
				return NONE;
		}

		result = numerator / denominator;
		return INTERSECT;
	}

	/**
	 * Calculates the intersection between 2 infinitely long, bidirectional
	 * lines.
	 */
	template<class T, class Scaler = typename T::Scalar>
	IntersectionType
	intersection(const T& p1, const T& v1, const T& p2, const T& v2, T& result,
		Scaler tolerance = std::numeric_limits<Scaler>::epsilon() * 4)
	{
		T      n2          = orthogonal(v2);
		Scaler denominator = v1.dot(n2);
		Scaler numerator   = (p2 - p1).dot(n2);

		if (std::abs(denominator) <= tolerance)
		{
			if (std::abs(numerator) <= tolerance)
				return COINCIDENT;
			else
				return NONE;
		}

		result = p1 + (v1 * (numerator / denominator));
		return INTERSECT;
	}

//	/**
//	 * Calculates the intersection between 2 infinitely long, bidirectional
//	 * lines, and the distance from p1 along the first line to the second line.
//	 * The distance is in units of ||v1||.
//	 * This is useful because the distance is calculated as an intermediary
//	 * anyway, so if you wanted this distance, and did not use this function you
//	 * would have to calculate the distance again.
//	 */
//	template<class T, int R, int O, int MR, int MC>
//	IntersectionType intersection(const Matrix<T, R, 1, O, MR, MC>& p1,
//		const Matrix<T, R, 1, O, MR, MC>& v1,
//		const Matrix<T, R, 1, O, MR, MC>& p2,
//		const Matrix<T, R, 1, O, MR, MC>& v2,
//		Matrix<T, R, 1, O, MR, MC>& result, T& distance,
//		T tolerance = std::numeric_limits<T>::epsilon() * 4)
//	{
//		Matrix<T, R, 1, O, MR, MC> n2          = orthogonal(v2);
//		T                          denominator = v1.dot(n2);
//		T                          numerator   = (p2 - p1).dot(n2);
//
//		if (std::abs(denominator) <= tolerance)
//		{
//			if (std::abs(numerator) <= tolerance)
//			{
//				distance = 0;
//				return COINCIDENT;
//			}
//			else
//				return NONE;
//		}
//
//		distance = numerator / denominator;
//		result   = p1 + v1 * distance;
//		return INTERSECT;
//	}

	/**
	 * Calculates the intersection between 2 infinitely long, bidirectional
	 * lines, and the distance from p1 along the first line to the second line.
	 * The distances are in units of ||v1||, and ||v2|| respectively.
	 * This is useful because the distance is calculated as an intermediary
	 * anyway, so if you wanted this distance, and did not use this function you
	 * would have to calculate the distance again.
	 */
	template<class T, class Scaler = typename T::Scalar>
	IntersectionType
	intersectionDistance(T& p1, T& v1, T& p2, T& v2, T& result, Scaler& d1,
		Scaler& d2,
		Scaler tolerance = std::numeric_limits<Scaler>::epsilon() * 4)
	{
		T      n2          = orthogonal(v2);
		Scaler denominator = v1.dot(n2);
		Scaler numerator   = (p2 - p1).dot(n2);

		if (std::abs(denominator) <= tolerance)
		{
			if (std::abs(numerator) <= tolerance)
			{
				d1 = d2 = 0;
				return COINCIDENT;
			}
			else
				return NONE;
		}

		d1     = numerator / denominator;
		result = p1 + v1 * d1;
		//TODO: No branch option
		for (uint32_t i = 0; i < T::RowsAtCompileTime; ++i)
			if (v2(i) != 0)
			{
				d2 = (result(i) - p2(i)) / v2(i);
				break;
			}
		return INTERSECT;
	}

//	/**
//	 * Calculates the intersection between the line defined by p and v, and the
//	 * axis aligned N-1 dimensional slice that is perpendicular to the given
//	 * dimension's axis.
//	 *
//	 * In 3D, this calculates the intersection between the line and either the
//	 * x, y, or z plane.
//	 * In 2D, this calculates the intersection between the line and either the x
//	 * or y axis.
//	 */
//	template<class T, int R, int O, int MR, int MC>
//	IntersectionType intersection(const Matrix<T, R, 1, O, MR, MC>& p,
//		const Matrix<T, R, 1, O, MR, MC>& v, const uint32_t dimension,
//		const T value, Matrix<T, R, 1, O, MR, MC>& result,
//		T tolerance = std::numeric_limits<T>::epsilon() * 4)
//	{
//		const T& denominator = v[dimension];
//		T numerator = value - p[dimension];
//
//		if (std::abs(denominator) <= tolerance)
//		{
//			if (std::abs(numerator) <= tolerance)
//				return COINCIDENT;
//			else
//				return NONE;
//		}
//
//		result = p + v * (numerator / denominator);
//		return INTERSECT;
//	}
//
//	/**
//	 * Calculates the intersection and distance between the line defined by p
//	 * and v, and the axis aligned N-1 dimensional slice that is perpendicular
//	 * to the given dimension's axis. The distance is in units of ||v||.
//	 *
//	 * In 3D, this calculates the intersection between the line and either the
//	 * x, y, or z plane.
//	 * In 2D, this calculates the intersection between the line and either the x
//	 * or y axis.
//	 */
//	template<class T, int R, int O, int MR, int MC>
//	IntersectionType intersection(const Matrix<T, R, 1, O, MR, MC>& p,
//		const Matrix<T, R, 1, O, MR, MC>& v, const uint32_t dimension,
//		const T value, Matrix<T, R, 1, O, MR, MC>& result, T& distance,
//		T tolerance = std::numeric_limits<T>::epsilon() * 4)
//	{
//		const T& denominator = v[dimension];
//		T numerator = value - p[dimension];
//
//		if (std::abs(denominator) <= tolerance)
//		{
//			if (std::abs(numerator) <= tolerance)
//			{
//				distance = 0;
//				return COINCIDENT;
//			}
//			else
//				return NONE;
//		}
//
//		distance = numerator / denominator;
//		result   = p + v * distance;
//		return INTERSECT;
//	}
//
//	/**
//	 * Calculates the distance between the line defined by p and v, and the axis
//	 * aligned N-1 dimensional slice that is perpendicular to the given
//	 * dimension's axis. The distance is in units of ||v||.
//	 *
//	 * In 3D, this calculates the distance between the line and either the x, y,
//	 * or z plane.
//	 * In 2D, this calculates the distance between the line and either the x or
//	 * y axis.
//	 */
//	template<class T, int R, int O, int MR, int MC>
//	IntersectionType distance(const Matrix<T, R, 1, O, MR, MC>& p,
//		const Matrix<T, R, 1, O, MR, MC>& v, const uint32_t dimension,
//		const T value, T& result,
//		T tolerance = std::numeric_limits<T>::epsilon() * 4)
//	{
//		const T& denominator = v[dimension];
//		T numerator = value - p[dimension];
//
//		if (std::abs(denominator) <= tolerance)
//		{
//			if (std::abs(numerator) <= tolerance)
//			{
//				result = 0;
//				return COINCIDENT;
//			}
//			else
//				return NONE;
//		}
//
//		result = numerator / denominator;
//		return INTERSECT;
//	}
}

#endif
