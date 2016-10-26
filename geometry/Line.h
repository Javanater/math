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
	class Line
	{
		private:
			typedef Line<DIM, ValueType> Lne;
			typedef Vector<DIM, ValueType> Vec;

		public:
			Vec point;
			Vec vector;

		public:
			Line() // Default
			{
			}

			Line(const Vec& point, const Vec& vector) :
					point(point), vector(vector) // Initializer
			{
			}

			Line(const Lne& line) :
					point(line.point), vector(line.vector) // Copy
			{
			}

			~Line()
			{
			}
	};

	typedef Line<2> Line2d;
	typedef Line<3> Line3d;
	typedef Line<4> Line4d;
}

#endif
