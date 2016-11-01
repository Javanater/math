/*
 * Line.h
 *
 *  Created on: Mar 13, 2016
 *      Author: Madison
 */

#ifndef LINE_H_
#define LINE_H_

#include "Vector.h"

namespace flabs
{
	template<uint32_t DIM, typename ValueType = double>
	class LineSegment
	{
		private:
			typedef LineSegment<DIM, ValueType>      Lne;
			typedef Eigen::Matrix<ValueType, DIM, 0> Vec;

		public:
			Vec start;
			Vec extends;

		public:
			LineSegment()
			{
			}

			LineSegment(std::initializer_list<const Vec> list)
			{
				auto it = list.begin();
				start = *it;
				extends = *(++it) - start;
			}

			LineSegment(const Vec start, const Vec finish) : start(start),
				extends(finish - start)
			{
			}

			LineSegment(const Lne& line) : start(line.start),
				extends(line.extends)
			{
			}

			~LineSegment()
			{
			}
	};

	typedef LineSegment<2> LineSegment2d;
	typedef LineSegment<3> LineSegment3d;
}

#endif
