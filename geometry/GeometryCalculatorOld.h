/*
 * GeometryCalculator.h
 *
 *  Created on: Mar 12, 2016
 *      Author: Madison
 */

#ifndef GEOMETRYCALCULATOR_H_
#define GEOMETRYCALCULATOR_H_

#include <limits>
#include "Vector.h"

namespace flabs
{
	enum IntersectionType
	{
		INTERSECT, COINCIDENT, NONE
	};

	std::ostream&
	operator<<(std::ostream& output, const IntersectionType& type);

	/**
	 * Calculates the intersection type between 2 infinitely long, bidirectional
	 * lines.
	 */
	template<uint32_t DIM, typename ValueType>
	IntersectionType intersection(const Vector<DIM, ValueType>& p1,
		const Vector<DIM, ValueType>& v1, const Vector<DIM, ValueType>& p2,
		const Vector<DIM, ValueType>& v2,
		ValueType tolerance = std::numeric_limits<ValueType>::epsilon() * 4)
	{
		Vector<DIM, ValueType> n2          = v2.orthogonal();
		ValueType              denominator = v1.dot(n2);

		if (std::abs(denominator) <= tolerance)
		{
			ValueType numerator = (p2 - p1).dot(n2);
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
	template<uint32_t DIM, typename ValueType>
	IntersectionType intersection(const Vector<DIM, ValueType>& p1,
		const Vector<DIM, ValueType>& v1, const Vector<DIM, ValueType>& p2,
		const Vector<DIM, ValueType>& v2, Vector<DIM, ValueType>& result,
		ValueType tolerance = std::numeric_limits<ValueType>::epsilon() * 4)
	{
		Vector<DIM, ValueType> n2          = v2.orthogonal();
		ValueType              denominator = v1.dot(n2);
		ValueType              numerator   = (p2 - p1).dot(n2);

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
	template<uint32_t DIM, typename ValueType>
	IntersectionType intersection(const Vector<DIM, ValueType>& p1,
		const Vector<DIM, ValueType>& v1, const Vector<DIM, ValueType>& p2,
		const Vector<DIM, ValueType>& v2, Vector<DIM, ValueType>& result,
		ValueType& distance,
		ValueType tolerance = std::numeric_limits<ValueType>::epsilon() * 4)
	{
		Vector<DIM, ValueType> n2          = v2.orthogonal();
		ValueType              denominator = v1.dot(n2);
		ValueType              numerator   = (p2 - p1).dot(n2);

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
	template<uint32_t DIM, typename ValueType>
	IntersectionType intersectionDistance(const Vector<DIM, ValueType>& p1,
		const Vector<DIM, ValueType>& v1, const Vector<DIM, ValueType>& p2,
		const Vector<DIM, ValueType>& v2, Vector<DIM, ValueType>& result,
		ValueType& d1, ValueType& d2,
		ValueType tolerance = std::numeric_limits<ValueType>::epsilon() * 4)
	{
		Vector<DIM, ValueType> n2          = v2.orthogonal();
		ValueType              denominator = v1.dot(n2);
		ValueType              numerator   = (p2 - p1).dot(n2);

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

		d1 = numerator / denominator;
		result   = p1 + v1 * d1;
		//TODO: No branch option
		for (uint32_t i = 0; i < DIM; ++i)
			if (v2[i] != 0)
			{
				d2 = (result[i] - p2[i]) / v2[i];
				break;
			}
		return INTERSECT;
	}

	/**
	 * Calculates the distance from p1 along the first line to the second line.
	 * This is useful for when you do not care about the intersection point,
	 * because it is never calculated. The distance is in units of ||v1||.
	 */
	template<uint32_t DIM, typename ValueType>
	IntersectionType
	distance(const Vector<DIM, ValueType>& p1, const Vector<DIM, ValueType>& v1,
		const Vector<DIM, ValueType>& p2, const Vector<DIM, ValueType>& v2,
		ValueType& result,
		ValueType tolerance = std::numeric_limits<ValueType>::epsilon() * 4)
	{
		Vector<DIM, ValueType> n2          = v2.orthogonal();
		ValueType              denominator = v1.dot(n2);
		ValueType              numerator   = (p2 - p1).dot(n2);

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
	template<uint32_t DIM, typename ValueType>
	IntersectionType intersection(const Vector<DIM, ValueType>& p,
		const Vector<DIM, ValueType>& v, const uint32_t dimension,
		const ValueType value, Vector<DIM, ValueType>& result,
		ValueType tolerance = std::numeric_limits<ValueType>::epsilon() * 4)
	{
		const ValueType& denominator = v[dimension];
		ValueType numerator = value - p[dimension];

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
	template<uint32_t DIM, typename ValueType>
	IntersectionType intersection(const Vector<DIM, ValueType>& p,
		const Vector<DIM, ValueType>& v, const uint32_t dimension,
		const ValueType value, Vector<DIM, ValueType>& result,
		ValueType& distance,
		ValueType tolerance = std::numeric_limits<ValueType>::epsilon() * 4)
	{
		const ValueType& denominator = v[dimension];
		ValueType numerator = value - p[dimension];

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
	template<uint32_t DIM, typename ValueType>
	IntersectionType
	distance(const Vector<DIM, ValueType>& p, const Vector<DIM, ValueType>& v,
		const uint32_t dimension, const ValueType value, ValueType& result,
		ValueType tolerance = std::numeric_limits<ValueType>::epsilon() * 4)
	{
		const ValueType& denominator = v[dimension];
		ValueType numerator = value - p[dimension];

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
