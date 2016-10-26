/*
 * SpatialTree.h
 *
 *  Created on: Mar 15, 2016
 *      Author: Madison
 */

#ifndef SPATIALTREE_H_
#define SPATIALTREE_H_

#include <stdlib.h>
#include <stdint.h>
#include "Vector.h"
#include "Ray.h"
#include "GeometryCalculator.h"

namespace flabs
{
	template<uint32_t DIM, class ExtendsVector, class VectorType = double>
	class SpatialTree
	{
		protected:
			typedef Vector<DIM, VectorType> Vec;
			typedef Ray<DIM, VectorType> Ry;
			typedef SpatialTree<DIM, ExtendsVector, VectorType> ST;
			ExtendsVector data;
			Vec corner;
			VectorType size;
			const uint32_t CHILD_COUNT = 1 << DIM;
			ST* children[CHILD_COUNT];

		public:
			SpatialTree(const ExtendsVector& data, const Vec& corner, const VectorType& size) :
					data(data), corner(corner), size(size)
			{
				for (uint32_t i = 0; i < CHILD_COUNT; ++i)
					children[i] = NULL;
			}

			virtual ~SpatialTree()
			{
				for (uint32_t i = 0; i < CHILD_COUNT; ++i)
					if (children[i] != NULL)
						delete children[i];
			}

			inline bool contains(const Ry& ray)
			{
				Vec t1 = (corner - ray.start) / ray.normalizedDirection;
				Vec t2 = (corner + size - ray.start) / ray.normalizedDirection;

				VectorType tmin = max(t1.min(), t2.min());
				VectorType tmax = min(t1.max(), t2.max());

				return tmax >= tmin;
			}

			inline bool bounds(const Vec& point)
			{
				return corner <= point && point < (corner + size);
			}

			inline bool insert(const ExtendsVector& data)
			{
				if (bounds(data))
				{
					insertNoCheck(data);
					return true;
				}
				else
					return false;
			}

		private:
			void insertNoCheck(const ExtendsVector& data)
			{
				const uint32_t index = getChildIndex(data);
				if (children[index] == NULL)
				{
					Vec corner(this->corner);
					for (int i = 0; i < DIM; ++i)
						if (data[i] - this->corner[i] >= size / 2)
							corner[i] += size / 2;
					children[index] = new ST(data, corner, size / 2);
				}
				else
					children[index]->insertNoCheck(data);
			}

			inline uint32_t getChildIndex(const Vec& point)
			{
				uint32_t index = 0;

				for (uint32_t i = 0; i < DIM; ++i)
				{
					index <<= 1;

					if (point.values[i] - corner.values[i] >= size / 2)
						index |= 1;
				}

				return index;
			}
	};
}

#endif
