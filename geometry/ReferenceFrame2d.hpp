//
// Created by Madison on 11/11/2016.
//

#ifndef PROJECTS_REFERENCEFRAME2D_HPP
#define PROJECTS_REFERENCEFRAME2D_HPP

#include <math.h>
#include "ReferenceFrame.hpp"

namespace flabs
{
	template<class ValueType = double>
	class ReferenceFrame2d : public ReferenceFrame<2, ValueType>
	{
		public:
			typedef typename ReferenceFrame<2, ValueType>::Tran Tran;

			ReferenceFrame2d(ValueType x = 0, ValueType y = 0, ValueType t = 0,
				ReferenceFrame2d* parent = nullptr)
			{
				this->parent = parent;

				setXOffset(x);
				setYOffset(y);
				setTOffset(t);

				this->transformationMatrix(2, 0) = 0;
				this->transformationMatrix(2, 1) = 0;
				this->transformationMatrix(2, 2) = 1;
			}

			inline ValueType getXOffset() const
			{
				return this->transformationMatrix(0, 2);
			}

			inline void setXOffset(ValueType x)
			{
				this->transformationMatrix(0, 2) = x;
			}

			inline ValueType getYOffset() const
			{
				return this->transformationMatrix(1, 2);
			}

			inline void setYOffset(ValueType y)
			{
				this->transformationMatrix(1, 2) = y;
			}

			inline ValueType getTOffset() const
			{
				return std::atan2(this->transformationMatrix(1, 0),
					this->transformationMatrix(1, 1));
			}

			inline void setTOffset(ValueType t)
			{
				this->transformationMatrix(0, 0) = std::cos(t);
				this->transformationMatrix(0, 1) = -std::sin(t);
				this->transformationMatrix(1, 0) = std::sin(t);
				this->transformationMatrix(1, 1) = std::cos(t);
			}

			inline void getXYT(ValueType& x, ValueType& y, ValueType& t) const
			{
				Tran offsetFromWorld = this->getOffsetFromWorld();
				x = offsetFromWorld(0, 2);
				y = offsetFromWorld(1, 2);
				t = std::atan2(offsetFromWorld(1, 0), offsetFromWorld(1, 1));
			}

			inline void getXY(ValueType& x, ValueType& y) const
			{
				Tran offsetFromWorld = this->getOffsetFromWorld();
				x = offsetFromWorld(0, 2);
				y = offsetFromWorld(1, 2);
			}
	};

	typedef ReferenceFrame2d<> ReferenceFrame2D;
}

#endif //PROJECTS_REFERENCEFRAME2D_HPP
