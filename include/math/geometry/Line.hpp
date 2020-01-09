/*
 * Line.h
 *
 *  Created on: Mar 13, 2016
 *      Author: Madison
 */

#ifndef LINE_H_
#define LINE_H_

#include "GeometryCalculator.hpp"

namespace flabs
{
	template<uint32_t DIM, typename ValueType = double>
	class Line
	{
		private:
			typedef Line<DIM, ValueType>             Lne;
			typedef Eigen::Matrix<ValueType, DIM, 0> Vec;

		public:
			Vec point;
			Vec vector;

		public:
			Line()
			{
			}

			Line(const Vec& point, const Vec& vector) :
				point(point), vector(vector)
			{
			}

			Line(const Lne& line) : point(line.point), vector(line.vector)
			{
			}

			~Line()
			{
			}

			inline IntersectionType intersect(const Lne& line, Vec& result)
			{
				return intersection(point, vector, line.point, line.vector,
					result);
			}

			inline ValueType distance(const Lne& line)
			{
				ValueType        distance;
				IntersectionType intersectionType =
									 distance(point, vector, line.point,
										 line.vector, distance);
				if (intersectionType == NONE)
					distance = std::numeric_limits<ValueType>::infinity();
				return distance;
			}
	};

	typedef Line<2> Line2d;
	typedef Line<3> Line3d;
	typedef Line<4> Line4d;
}

#endif
