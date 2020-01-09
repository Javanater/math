//
// Created by Madison on 8/31/2016.
//

#ifndef PROJECTS_POSE_HPP
#define PROJECTS_POSE_HPP

#include <Eigen/Eigen>
#include "ReferenceFrame.hpp"

namespace flabs
{
template<uint32_t DIM, class ValueType = double>
class ReferencePoint
{
public:
	typedef ReferenceFrame <DIM, ValueType>            Ref;
	typedef Eigen::Matrix<ValueType, DIM, 1>           Vec;
	typedef Eigen::Matrix<ValueType, DIM, DIM>         Rot;
	typedef Eigen::Matrix<ValueType, DIM + 1, DIM + 1> Tran;

public:
	Ref* parent;
	Tran transformationMatrix;

public:
	ReferencePoint(Ref* parent = nullptr) : parent(parent)
	{
		transformationMatrix.setIdentity();
	}

	template<uint32_t U = DIM, class = typename std::enable_if<DIM == 2>::type>
	ReferencePoint(ValueType x, ValueType y, ValueType yaw,
		Ref* parent = nullptr) : parent(parent)
	{
		setXOffset(x);
		setYOffset(y);
		setYawOffset(yaw);

		transformationMatrix(2, 0) = 0;
		transformationMatrix(2, 1) = 0;
		transformationMatrix(2, 2) = 1;
	}

	virtual ~ReferencePoint()
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

	Vec getTranslationOffset()
	{
		return transformationMatrix.block<DIM, 1>(0, 2);
	}

	Vec getWorldPosition()
	{
		if (parent)
			return (parent->getOffsetFromWorld() * transformationMatrix)
				.block<DIM, 1>(0, 2);
		else
			return transformationMatrix.block<DIM, 1>(0, 2);
	}

	inline ValueType getXOffset() const
	{
		return transformationMatrix(0, DIM);
	};

	inline Ref& setXOffset(ValueType x)
	{
		transformationMatrix(0, DIM) = x;
		return *this;
	};

	template<uint32_t Dummy = DIM>
	inline typename std::enable_if<Dummy >= 2, ValueType>::type
	getYOffset() const
	{
		return transformationMatrix(1, DIM);
	};

	template<uint32_t Dummy = DIM>
	inline typename std::enable_if<Dummy >= 2, Ref&>::type
	setYOffset(ValueType y)
	{
		transformationMatrix(1, DIM) = y;
		return *this;
	};

	template<uint32_t Dummy = DIM>
	inline typename std::enable_if<Dummy >= 3, ValueType>::type
	getZOffset() const
	{
		return transformationMatrix(2, DIM);
	};

	template<uint32_t Dummy = DIM>
	inline typename std::enable_if<Dummy >= 3, Ref&>::type
	setZOffset(ValueType z)
	{
		transformationMatrix(2, DIM) = z;
		return *this;
	};

	template<uint32_t Dummy = DIM>
	inline typename std::enable_if<Dummy == 2, ValueType>::type
	getYawOffset() const
	{
		return std::atan2(this->transformationMatrix(1, 0),
			this->transformationMatrix(1, 1));
	};

	template<uint32_t Dummy = DIM>
	inline typename std::enable_if<Dummy == 2, Ref&>::type
	setYawOffset(ValueType yaw)
	{
		transformationMatrix(0, 0) = std::cos(yaw);
		transformationMatrix(0, 1) = -std::sin(yaw);
		transformationMatrix(1, 0) = std::sin(yaw);
		transformationMatrix(1, 1) = std::cos(yaw);
		return *this;
	};

	template<uint32_t Dummy = DIM>
	inline typename std::enable_if<Dummy == 2, const Ref&>::type
	getXYYaw(ValueType& x, ValueType& y, ValueType& yaw) const
	{
		Tran offsetFromWorld = getOffsetFromWorld();
		x   = offsetFromWorld(0, 2);
		y   = offsetFromWorld(1, 2);
		yaw = std::atan2(offsetFromWorld(1, 0), offsetFromWorld(1, 1));
		return *this;
	};

	template<uint32_t Dummy = DIM>
	inline typename std::enable_if<Dummy == 2, Ref&>::type
	getXYYaw(ValueType& x, ValueType& y, ValueType& yaw)
	{
		Tran offsetFromWorld = getOffsetFromWorld();
		x   = offsetFromWorld(0, 2);
		y   = offsetFromWorld(1, 2);
		yaw = std::atan2(offsetFromWorld(1, 0), offsetFromWorld(1, 1));
		return *this;
	};

	template<uint32_t Dummy = DIM>
	inline typename std::enable_if<Dummy == 2, const Ref&>::type
	getXY(ValueType& x, ValueType& y) const
	{
		Tran offsetFromWorld = getOffsetFromWorld();
		x = offsetFromWorld(0, 2);
		y = offsetFromWorld(1, 2);
		return *this;
	};

	template<uint32_t Dummy = DIM>
	inline typename std::enable_if<Dummy == 2, Ref&>::type
	getXY(ValueType& x, ValueType& y)
	{
		Tran offsetFromWorld = getOffsetFromWorld();
		x = offsetFromWorld(0, 2);
		y = offsetFromWorld(1, 2);
		return *this;
	};
};

typedef ReferenceFrame<2, double> ReferenceFrame2d;
typedef ReferenceFrame<3, double> ReferenceFrame3d;
}

#endif //PROJECTS_POSE_HPP
