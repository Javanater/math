//
// Created by Madison on 8/31/2016.
//

#ifndef PROJECTS_POSE_HPP
#define PROJECTS_POSE_HPP

#include <cstdint>
#include "../matrix/Matrix.h"

namespace flabs
{
	template<uint32_t DIM, class TYPE>
	class Pose
	{
		private:
			Matrix<DIM * 2 + 1, DIM * 2 + 1> poseMatrix;

		public:
			Pose() : poseMatrix(0)
			{
			}
	};
}

#endif //PROJECTS_POSE_HPP
