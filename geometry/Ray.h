/*
 * Line.h
 *
 *  Created on: Mar 13, 2016
 *      Author: Madison
 */

#ifndef RAY_H_
#define RAY_H_

#include "GeometryCalculator.h"
#include "LineSegment.h"

namespace flabs
{
	template<uint32_t DIM, typename ValueType = double>
	class Ray
	{
		private:
			typedef Ray<DIM, ValueType>              Ry;
			typedef LineSegment<DIM, ValueType>      Seg;
			typedef Eigen::Matrix<ValueType, DIM, 0> Vec;

		public:
			Vec start;
			Vec normalizedDirection;

		public:
			Ray()
			{
			}

			Ray(const Vec start, const Vec normalizedDirection) :
				start(start), normalizedDirection(normalizedDirection)
			{
			}

			Ray(const Ry& ray) :
				start(ray.start), normalizedDirection(ray.normalizedDirection)
			{
			}

			~Ray()
			{
			}

			IntersectionType
			distance(const LineSegment<DIM, ValueType>& segment,
				ValueType& distance) const
			{
				Vec              intersection;
				ValueType        d2;
				IntersectionType intersectionType = intersectionDistance(start,
					normalizedDirection, segment.start, segment.extends,
					intersection, distance, d2);
				if (d2 < 0 || d2 > 1 || distance < 0)
					return NONE;
				else
					return intersectionType;
			}
	};

	typedef Ray<2> Ray2d;
	typedef Ray<3> Ray3d;
}

#endif
