/*
 * GeometryCalculator.cpp
 *
 *  Created on: Apr 19, 2016
 *      Author: Madison
 */

#include "GeometryCalculator.h"

namespace flabs
{
	std::ostream& operator<<(std::ostream& output, const IntersectionType& type)
	{
		switch (type)
		{
			case INTERSECT:
				output << "INTERSECT";
				break;
			case COINCIDENT:
				output << "COINCIDENT";
				break;
			case NONE:
				output << "NONE";
				break;
			default:
				output << "UNKNOWN_IntersectionType";
				break;
		}
		return output;
	}
}
