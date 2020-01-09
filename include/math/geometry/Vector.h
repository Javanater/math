/*
 * Vector.h
 *
 *  Created on: Mar 11, 2016
 *      Author: Madison
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include <initializer_list>
#include <cmath>
#include <iostream>
#include <functional>
#include <limits>

namespace flabs
{
	template<uint32_t DIM, typename ValueType = double>
	class Vector
	{
		private:
			typedef Vector<DIM, ValueType> Vec;
			ValueType                      values[DIM];

		public:
			Vector() // Default
			{
			}

			Vector(const Vec& vector) // Copy
			{
				for (uint32_t i = 0; i < DIM; ++i)
					values[i] = vector.values[i];
			}

			Vector(ValueType value) // Initializer
			{
				for (uint32_t i = 0; i < DIM; ++i)
					values[i] = value;
			}

			Vector(const ValueType values[DIM]) // Initializer
			{
				for (uint32_t i = 0; i < DIM; ++i)
					this->values[i] = values[i];
			}

			Vector(std::initializer_list<ValueType> list) // Initializer
			{
				auto          it = list.begin();
				for (uint32_t i  = 0; i < DIM; ++i)
					values[i] = *it++;
			}

			Vector(std::function<ValueType()>& generator)
			{
				for (uint32_t i = 0; i < DIM; ++i)
					values[i] = generator();
			}

			static Vec fromGenerator(std::function<ValueType()>& generator)
			{
				Vec           v;
				for (uint32_t i = 0; i < DIM; ++i)
					v.values[i] = generator();
				return v;
			}

			inline Vec orthogonal() const
			{
				if (DIM == 2)
					return orthogonal2d();
				else if(DIM == 3)
					return orthogonal3d();
				else
				{
					Vec           normal;
					for (uint32_t skip = 0; skip < DIM; ++skip)
					{
						ValueType     sum = 0;
						for (uint32_t i   = 0; i < DIM; ++i)
						{
							if (i != skip)
							{
								normal[i] = values[skip];
								sum -= values[i];
							}
						}
						normal[skip] = sum;
						if (sum != 0 && values[skip] != 0)
							break;
					}
					return normal;
				}
			}

		private:
			//TODO: Templatize, so all dimensions are optimized at compile time.
			inline Vec orthogonal2d() const
			{
				Vec           normal;
				normal[0] = -values[1];
				normal[1] = values[0];
				return normal;
			}

			inline Vec orthogonal3d() const
			{
				Vec normal[2];
				normal[0][0] = -values[1] - values[2];
				normal[0][1] = values[0];
				normal[0][2] = values[0];
				normal[1][0] = -values[0] - values[2];
				normal[1][1] = values[1];
				normal[1][2] = values[1];
				int select = values[1] == 0 && values[0] == -values[2];
				return normal[select];
			}

		public:
			~Vector()
			{
			}

			friend std::ostream&
			operator<<(std::ostream& output, const Vec& vector)
			{
				output << "<";
				if (DIM > 0)
				{
					output << vector.values[0];
					for (uint32_t i = 1; i < DIM; ++i)
						output << "," << vector.values[i];
				}
				output << ">";
				return output;
			}

			inline Vec& operator=(const Vec& vector)
			{
				for (uint32_t i = 0; i < DIM; ++i)
					values[i] = vector.values[i];
				return *this;
			}

			inline Vec& operator=(std::initializer_list<ValueType> list)
			{
				auto          it = list.begin();
				for (uint32_t i  = 0; i < DIM; ++i)
					values[i] = *it++;
				return *this;
			}

			inline ValueType dot(const Vec& vector) const
			{
				ValueType     dotProduct = 0;
				for (uint32_t i          = 0; i < DIM; ++i)
					dotProduct += values[i] * vector.values[i];
				return dotProduct;
			}

			inline ValueType norm() const
			{
				ValueType     norm = 0;
				for (uint32_t i    = 0; i < DIM; ++i)
					norm += values[i] * values[i];
				return std::sqrt(norm);
			}

			inline Vec& normalize()
			{
				return operator/=(norm());
			}

			inline Vec& operator+=(const Vec& vector)
			{
				for (uint32_t i = 0; i < DIM; ++i)
					values[i] += vector.values[i];
				return *this;
			}

			inline Vec& operator-=(const Vec& vector)
			{
				for (uint32_t i = 0; i < DIM; ++i)
					values[i] -= vector.values[i];
				return *this;
			}

			inline Vec operator*=(const Vec& vector)
			{
				for (uint32_t i = 0; i < DIM; ++i)
					values[i] *= vector.values[i];
				return *this;
			}

			inline Vec operator/=(const Vec& vector)
			{
				for (uint32_t i = 0; i < DIM; ++i)
					values[i] /= vector.values[i];
				return *this;
			}

			inline Vec operator+(Vec vector) const
			{
				for (uint32_t i = 0; i < DIM; ++i)
					vector.values[i] += values[i];
				return vector;
			}

			inline Vec operator-(Vec vector) const
			{
				for (uint32_t i = 0; i < DIM; ++i)
					vector.values[i] = values[i] - vector.values[i];
				return vector;
			}

			inline Vec operator*(Vec vector) const
			{
				for (uint32_t i = 0; i < DIM; ++i)
					vector.values[i] = values[i] * vector.values[i];
				return vector;
			}

			inline Vec operator/(Vec vector) const
			{
				for (uint32_t i = 0; i < DIM; ++i)
					vector.values[i] = values[i] / vector.values[i];
				return vector;
			}

			inline bool operator==(const Vec& vector) const
			{
				for (uint32_t i = 0; i < DIM; ++i)
					if (values[i] != vector.values[i])
						return false;
				return true;
			}

			inline bool equals(const Vec& vector,
				double percentageTolerance = .00000001) const
			{
				for (uint32_t i = 0; i < DIM; ++i)
					if (std::abs(values[i] - vector.values[i]) / values[i] >
						percentageTolerance)
						return false;
				return true;
			}

			inline Vec& operator+=(const ValueType scaler)
			{
				for (uint32_t i = 0; i < DIM; ++i)
					values[i] += scaler;
				return *this;
			}

			inline Vec& operator-=(const ValueType scaler)
			{
				for (uint32_t i = 0; i < DIM; ++i)
					values[i] -= scaler;
				return *this;
			}

			inline Vec& operator*=(const ValueType scaler)
			{
				for (uint32_t i = 0; i < DIM; ++i)
					values[i] *= scaler;
				return *this;
			}

			inline Vec& operator/=(const ValueType scaler)
			{
				for (uint32_t i = 0; i < DIM; ++i)
					values[i] /= scaler;
				return *this;
			}

			inline ValueType& operator[](uint32_t index)
			{
				return values[index];
			}

			inline const ValueType& operator[](uint32_t index) const
			{
				return values[index];
			}

			friend inline Vec operator+(const ValueType scaler, Vec vector)
			{
				return vector + scaler;
			}

			friend inline Vec operator+(Vec vector, const ValueType scaler)
			{
				for (uint32_t i = 0; i < DIM; ++i)
					vector.values[i] += scaler;
				return vector;
			}

			inline Vec operator-() const
			{
				Vec           vector;
				for (uint32_t i = 0; i < DIM; ++i)
					vector.values[i] = -values[i];
				return vector;
			}

			friend inline Vec operator-(Vec vector, const ValueType scaler)
			{
				for (uint32_t i = 0; i < DIM; ++i)
					vector.values[i] -= scaler;
				return vector;
			}

			friend inline Vec operator*(const ValueType scaler, Vec vector)
			{
				return vector * scaler;
			}

			friend inline Vec operator*(Vec vector, const ValueType scaler)
			{
				for (uint32_t i = 0; i < DIM; ++i)
					vector.values[i] *= scaler;
				return vector;
			}

			friend inline Vec operator/(Vec vector, const ValueType scaler)
			{
				for (uint32_t i = 0; i < DIM; ++i)
					vector.values[i] /= scaler;
				return vector;
			}

			inline ValueType min()
			{
				ValueType     min = values[0];
				for (uint32_t i   = 1; i < DIM; ++i)
					min = min(min, values[i]);
				return min;
			}

			inline ValueType max()
			{
				ValueType     max = values[0];
				for (uint32_t i   = 1; i < DIM; ++i)
					max = max(max, values[i]);
				return max;
			}

			inline bool operator<(const Vec& vector) const
			{
				for (uint32_t i = 1; i < DIM; ++i)
					if (values[i] >= vector.values[i])
						return false;
				return true;
			}

			inline bool operator<=(const Vec& vector) const
			{
				for (uint32_t i = 1; i < DIM; ++i)
					if (values[i] > vector.values[i])
						return false;
				return true;
			}

			inline bool operator>(const Vec& vector) const
			{
				for (uint32_t i = 1; i < DIM; ++i)
					if (values[i] <= vector.values[i])
						return false;
				return true;
			}

			inline bool operator>=(const Vec& vector) const
			{
				for (uint32_t i = 1; i < DIM; ++i)
					if (values[i] < vector.values[i])
						return false;
				return true;
			}

			inline ValueType* ptr()
			{
				return (ValueType*) values;
			}

			inline const ValueType* ptr() const
			{
				return (ValueType*) values;
			}
	};

	//TODO: Specialize std::abs. Not possible???
//	namespace std
//	{
//		template<uint32_t DIM, typename ValueType>
//		template<>
//		Vector<DIM, ValueType> abs<Vector<DIM, ValueType>>(
//			const Vector<DIM, ValueType>& vector)
//		{
//			Vector<DIM, ValueType> absVector;
//			for (uint32_t i = 0; i < DIM; ++i)
//				absVector[i] = abs(absVector[i]);
//			return absVector;
//		}
//	}

	typedef Vector<2> Vector2d;
	typedef Vector<3> Vector3d;
	typedef Vector<4> Vector4d;
}

#endif
