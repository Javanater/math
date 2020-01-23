/*
 * Vector.h
 *
 *  Created on: Mar 11, 2016
 *      Author: Madison
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include <Eigen/Eigen>

namespace flabs
{
template<int DIM, class ValueType = double>
using Vector = Eigen::Matrix<ValueType, DIM, 1>;

using Vector2d = Vector<2, double>;
using Vector3d = Vector<3, double>;
using Vector4d = Vector<4, double>;
}

#endif
