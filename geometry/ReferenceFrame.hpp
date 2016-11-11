//
// Created by Madison on 8/31/2016.
//

#ifndef PROJECTS_POSE_HPP
#define PROJECTS_POSE_HPP

#include <Eigen/Eigen>

namespace flabs
{
	template<uint32_t DIM, class ValueType = double>
	class ReferenceFrame
	{
		public:
			typedef ReferenceFrame<DIM, ValueType>             Ref;
			typedef Eigen::Matrix<ValueType, DIM, 1>           Vec;
			typedef Eigen::Matrix<ValueType, DIM, DIM>         Rot;
			typedef Eigen::Matrix<ValueType, DIM + 1, DIM + 1> Tran;

		protected:
			Ref* parent;
			Tran transformationMatrix;

		public:
			ReferenceFrame() : parent(nullptr)
			{
			}

			Tran getOffsetFromWorld() const
			{
				if (parent)
					return parent->getOffsetFromWorld() * transformationMatrix;
				else
					return transformationMatrix;
			}

			Ref* getParent() const
			{
				return parent;
			}

			void setParent(Ref* parent)
			{
				ReferenceFrame::parent = parent;
			}
	};
}

#endif //PROJECTS_POSE_HPP
