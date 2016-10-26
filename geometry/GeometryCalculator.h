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

namespace flabs
{
	enum IntersectionType
	{
		INTERSECT, COINCIDENT, NONE
	};

	std::ostream&
	operator<<(std::ostream& output, const IntersectionType& type);

	template<class T>
	inline Eigen::Matrix<T, Eigen::MatrixBase<T>::RowsAtCompileTime, 1>
	orthogonal2d(Eigen::MatrixBase<T> v)
	{
		std::swap(v(0), v(1));
		v(0) = -v(0);
		return v;
	}

	template<class T>
	inline Eigen::Matrix<T, Eigen::MatrixBase<T>::RowsAtCompileTime, 1>
	orthogonal3d(const Eigen::MatrixBase<T>& v)
	{
		Eigen::Matrix<T, Eigen::MatrixBase<T>::RowsAtCompileTime, 1> normal[2];
		normal[0](0) = -v(1) - v(2);
		normal[0](1) = v(0);
		normal[0](2) = v(0);
		normal[1](0) = -v(0) - v(2);
		normal[1](1) = v(1);
		normal[1](2) = v(1);
		int select = v(1) == 0 && v(0) == -v(2);
		return normal[select];
	}

	//TODO: Templatize, so all dimensions are optimized at compile time.
	template<class T>
	inline Eigen::Matrix<T, Eigen::MatrixBase<T>::RowsAtCompileTime, 1>
	orthogonal(const Eigen::MatrixBase<T>& v)
	{
		if (Eigen::MatrixBase<T>::RowsAtCompileTime == 2)
			return orthogonal2d(v);
		else if (Eigen::MatrixBase<T>::RowsAtCompileTime == 3)
			return orthogonal3d(v);
		else
		{
			Eigen::Matrix<T, Eigen::MatrixBase<T>::RowsAtCompileTime, 1> normal;

			for (int skip = 0; skip < Eigen::MatrixBase<T>::RowsAtCompileTime;
				++skip)
			{
				T        sum = 0;
				for (int i   = 0; i < Eigen::MatrixBase<T>::RowsAtCompileTime;
					++i)
				{
					if (i != skip)
					{
						normal(i) = v(skip);
						sum -= v(i);
					}
				}
				normal(skip) = sum;
				if (sum != 0 && v(skip) != 0)
					break;
			}
			return normal;
		}
	}

	/**
	 * Calculates the intersection type between 2 infinitely long, bidirectional
	 * lines.
	 */
	template<class T1, class T2, class T3, class T4>
	IntersectionType intersection(const Eigen::MatrixBase<T1>& p1,
		const Eigen::MatrixBase<T2>& v1, const Eigen::MatrixBase<T3>& p2,
		const Eigen::MatrixBase<T4>& v2,
		T1::Scaler tolerance = std::numeric_limits<T1::Scaler>::epsilon() * 4)
	{
		Eigen::Matrix<T2, Eigen::MatrixBase<T2>::RowsAtCompileTime, 1>
				   n2          = orthogonal(v2);
		T1::Scaler denominator = v1.dot(n2);

		if (std::abs(denominator) <= tolerance)
		{
			T1::Scaler numerator = (p2 - p1).dot(n2);
			if (std::abs(numerator) <= tolerance)
				return COINCIDENT;
			else
				return NONE;
		}

		return INTERSECT;
	}

	/**
	 * Calculates the intersection between 2 infinitely long, bidirectional
	 * lines.
	 */
	template<class T, int R, int O, int MR, int MC>
	IntersectionType intersection(const Matrix<T, R, 1, O, MR, MC>& p1,
		const Matrix<T, R, 1, O, MR, MC>& v1,
		const Matrix<T, R, 1, O, MR, MC>& p2,
		const Matrix<T, R, 1, O, MR, MC>& v2,
		Matrix<T, R, 1, O, MR, MC>& result,
		T tolerance = std::numeric_limits<T>::epsilon() * 4)
	{
		Matrix<T, R, 1, O, MR, MC> n2          = orthogonal(v2);
		T                          denominator = v1.dot(n2);
		T                          numerator   = (p2 - p1).dot(n2);

		if (std::abs(denominator) <= tolerance)
		{
			if (std::abs(numerator) <= tolerance)
				return COINCIDENT;
			else
				return NONE;
		}

		result = p1 + v1 * (numerator / denominator);
		return INTERSECT;
	}

	/**
	 * Calculates the intersection between 2 infinitely long, bidirectional
	 * lines, and the distance from p1 along the first line to the second line.
	 * The distance is in units of ||v1||.
	 * This is useful because the distance is calculated as an intermediary
	 * anyway, so if you wanted this distance, and did not use this function you
	 * would have to calculate the distance again.
	 */
	template<class T, int R, int O, int MR, int MC>
	IntersectionType intersection(const Matrix<T, R, 1, O, MR, MC>& p1,
		const Matrix<T, R, 1, O, MR, MC>& v1,
		const Matrix<T, R, 1, O, MR, MC>& p2,
		const Matrix<T, R, 1, O, MR, MC>& v2,
		Matrix<T, R, 1, O, MR, MC>& result, T& distance,
		T tolerance = std::numeric_limits<T>::epsilon() * 4)
	{
		Matrix<T, R, 1, O, MR, MC> n2          = orthogonal(v2);
		T                          denominator = v1.dot(n2);
		T                          numerator   = (p2 - p1).dot(n2);

		if (std::abs(denominator) <= tolerance)
		{
			if (std::abs(numerator) <= tolerance)
			{
				distance = 0;
				return COINCIDENT;
			}
			else
				return NONE;
		}

		distance = numerator / denominator;
		result   = p1 + v1 * distance;
		return INTERSECT;
	}

	/**
	 * Calculates the intersection between 2 infinitely long, bidirectional
	 * lines, and the distance from p1 along the first line to the second line.
	 * The distances are in units of ||v1||, and ||v2|| respectively.
	 * This is useful because the distance is calculated as an intermediary
	 * anyway, so if you wanted this distance, and did not use this function you
	 * would have to calculate the distance again.
	 */
	template<class T, int R, int O, int MR, int MC>
	IntersectionType intersectionDistance(const Matrix<T, R, 1, O, MR, MC>& p1,
		const Matrix<T, R, 1, O, MR, MC>& v1,
		const Matrix<T, R, 1, O, MR, MC>& p2,
		const Matrix<T, R, 1, O, MR, MC>& v2,
		Matrix<T, R, 1, O, MR, MC>& result, T& d1, T& d2,
		T tolerance = std::numeric_limits<T>::epsilon() * 4)
	{
		Matrix<T, R, 1, O, MR, MC> n2          = orthogonal(v2);
		T                          denominator = v1.dot(n2);
		T                          numerator   = (p2 - p1).dot(n2);

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
		for (uint32_t i = 0; i < R; ++i)
			if (v2(i) != 0)
			{
				d2 = (result(i) - p2(i)) / v2(i);
				break;
			}
		return INTERSECT;
	}

	/**
	 * Calculates the distance from p1 along the first line to the second line.
	 * This is useful for when you do not care about the intersection point,
	 * because it is never calculated. The distance is in units of ||v1||.
	 */
	template<class T, int R, int O, int MR, int MC>
	IntersectionType distance(const Matrix<T, R, 1, O, MR, MC>& p1,
		const Matrix<T, R, 1, O, MR, MC>& v1,
		const Matrix<T, R, 1, O, MR, MC>& p2,
		const Matrix<T, R, 1, O, MR, MC>& v2, T& result,
		T tolerance = std::numeric_limits<T>::epsilon() * 4)
	{
		Matrix<T, R, 1, O, MR, MC> n2          = orthogonal(v2);
		T                          denominator = v1.dot(n2);
		T                          numerator   = (p2 - p1).dot(n2);

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
	 * Calculates the intersection between the line defined by p and v, and the
	 * axis aligned N-1 dimensional slice that is perpendicular to the given
	 * dimension's axis.
	 *
	 * In 3D, this calculates the intersection between the line and either the
	 * x, y, or z plane.
	 * In 2D, this calculates the intersection between the line and either the x
	 * or y axis.
	 */
	template<class T, int R, int O, int MR, int MC>
	IntersectionType intersection(const Matrix<T, R, 1, O, MR, MC>& p,
		const Matrix<T, R, 1, O, MR, MC>& v, const uint32_t dimension,
		const T value, Matrix<T, R, 1, O, MR, MC>& result,
		T tolerance = std::numeric_limits<T>::epsilon() * 4)
	{
		const T& denominator = v[dimension];
		T numerator = value - p[dimension];

		if (std::abs(denominator) <= tolerance)
		{
			if (std::abs(numerator) <= tolerance)
				return COINCIDENT;
			else
				return NONE;
		}

		result = p + v * (numerator / denominator);
		return INTERSECT;
	}

	/**
	 * Calculates the intersection and distance between the line defined by p
	 * and v, and the axis aligned N-1 dimensional slice that is perpendicular
	 * to the given dimension's axis. The distance is in units of ||v||.
	 *
	 * In 3D, this calculates the intersection between the line and either the
	 * x, y, or z plane.
	 * In 2D, this calculates the intersection between the line and either the x
	 * or y axis.
	 */
	template<class T, int R, int O, int MR, int MC>
	IntersectionType intersection(const Matrix<T, R, 1, O, MR, MC>& p,
		const Matrix<T, R, 1, O, MR, MC>& v, const uint32_t dimension,
		const T value, Matrix<T, R, 1, O, MR, MC>& result, T& distance,
		T tolerance = std::numeric_limits<T>::epsilon() * 4)
	{
		const T& denominator = v[dimension];
		T numerator = value - p[dimension];

		if (std::abs(denominator) <= tolerance)
		{
			if (std::abs(numerator) <= tolerance)
			{
				distance = 0;
				return COINCIDENT;
			}
			else
				return NONE;
		}

		distance = numerator / denominator;
		result   = p + v * distance;
		return INTERSECT;
	}

	/**
	 * Calculates the distance between the line defined by p and v, and the axis
	 * aligned N-1 dimensional slice that is perpendicular to the given
	 * dimension's axis. The distance is in units of ||v||.
	 *
	 * In 3D, this calculates the distance between the line and either the x, y,
	 * or z plane.
	 * In 2D, this calculates the distance between the line and either the x or
	 * y axis.
	 */
	template<class T, int R, int O, int MR, int MC>
	IntersectionType distance(const Matrix<T, R, 1, O, MR, MC>& p,
		const Matrix<T, R, 1, O, MR, MC>& v, const uint32_t dimension,
		const T value, T& result,
		T tolerance = std::numeric_limits<T>::epsilon() * 4)
	{
		const T& denominator = v[dimension];
		T numerator = value - p[dimension];

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
}

#endif
