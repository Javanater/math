/*******************************************************************************
 File:				LUDecomposition.h

 Author: 			Gaspard Petit (gaspardpetit@gmail.com)
 Last Revision:		March 14, 2007

 This code may be reused without my permission as long as credits are given to
 the original author.  If you find bugs, please send me a note...
*******************************************************************************/
#ifndef __MK_GEOMETRY_LUDECOMPOSITION__
#define __MK_GEOMETRY_LUDECOMPOSITION__

#include "../geometry/Vector.h"

//==============================================================================
// Factorise a matrix A into an L(ower) and U(pper) matrix such that
// A = L*U or P*A = L*U where P is a pivot matrix (changes the row order).
//==============================================================================

namespace flabs
{
	template<uint32_t ROWS, uint32_t COLS, class TYPE>
	class Matrix;

	class LUDecomposition
	{
		public:
			template<typename _T, uint32_t _SIZE>
			static bool decompose(const Matrix<_SIZE, _SIZE, _T>& matrix,
				Matrix<_SIZE, _SIZE, _T>* l, Matrix<_SIZE, _SIZE, _T>* u)
			{
				*u = matrix;
				*l = Matrix<_SIZE, _SIZE, _T>::identity();
				return LUSubDecomposition<_T, 0, _SIZE, _SIZE>::decomposeRecur(
					u->ptr(), l->ptr());
			}

			template<typename _T, uint32_t _SIZE>
			static bool decompose(const Matrix<_SIZE, _SIZE, _T>& matrix,
				Matrix<_SIZE, _SIZE, _T>* l, Matrix<_SIZE, _SIZE, _T>* u,
				Vector<_SIZE, int>* p)
			{
				for (uint32_t i = 0; i < _SIZE; ++i)
					(*p)[i] = i;
				*u = matrix;
				*l = Matrix<_SIZE, _SIZE, _T>::identity();
				return LUSubDecomposition<_T, 0, _SIZE, _SIZE>::decomposeRecur(
					u->ptr(), l->ptr(), p->ptr());
			}

			template<typename _T, uint32_t _SIZE, uint32_t _BXWIDTH>
			static bool solve(const Matrix<_SIZE, _SIZE, _T>& l,
				const Matrix<_SIZE, _SIZE, _T>& u,
				Matrix<_SIZE, _BXWIDTH, _T>* xb)
			{
				return LUSubDecomposition<_T, 0, _SIZE, _SIZE>::solve(l, u, xb);
			}

		private:
			template<
				typename _T, uint32_t _OFFSET, uint32_t _WIDTH, uint32_t _HEIGHT>
			class LUSubDecomposition
			{
				public:
					static bool decomposeRecur(_T* u, _T* l)
					{
						if (!decompose(u, l))
							return false;

						return LUSubDecomposition<_T, _OFFSET + 1, _WIDTH,
							_HEIGHT>::decomposeRecur(u, l);
					}

					static bool decompose(_T* u, _T* l)
					{
						const uint32_t width  = _WIDTH;
						const uint32_t height = _HEIGHT;

						// normalize the row
						_T factor = u[_OFFSET * width + _OFFSET];
						l[_OFFSET * width + _OFFSET] = factor;

						u[_OFFSET * width + _OFFSET] = 1;
						RowOperation<_T, _WIDTH - _OFFSET - 1>::multiply(
							&u[_OFFSET * width + _OFFSET + 1],
							&u[_OFFSET * width + _OFFSET + 1],
							_T(1.0) / factor);

						for (uint32_t j = _OFFSET + 1; j < height; ++j)
						{
							factor = u[j * width + _OFFSET];
							l[j * width + _OFFSET] = factor;
							u[j * width + _OFFSET] = 0;
							RowOperation<_T,
								_WIDTH - _OFFSET - 1>::addRowTimesFactor(
								&u[j * width + _OFFSET + 1],
								&u[j * width + _OFFSET + 1],
								&u[_OFFSET * width + _OFFSET + 1], -factor);
						}

						return true;
					}

					static bool decomposeRecur(_T* u, _T* l, int* p)
					{
						if (!decompose(u, l, p))
							return false;

						return LUSubDecomposition<_T, _OFFSET + 1, _WIDTH,
							_HEIGHT>::decomposeRecur(u, l, p);
					}

					static bool decompose(_T* u, _T* l, int* p)
					{
						const uint32_t width  = _WIDTH;
						const uint32_t height = _HEIGHT;

						//	swap with a lower row so that we have the highest value factor
						_T            max    = std::abs(
							u[_OFFSET * width + _OFFSET]);
						uint32_t      maxRow = _OFFSET;
						for (uint32_t j      = _OFFSET + 1; j < height; ++j)
						{
							_T v = std::abs(u[j * width + _OFFSET]);
							if (v > max)
							{
								max    = v;
								maxRow = j;
							}
						}

						if (maxRow != _OFFSET)
						{
							int temp = p[_OFFSET];
							p[_OFFSET] = p[maxRow];
							p[maxRow]  = temp;
							RowOperation<_T, _WIDTH - _OFFSET>::swap(
								&u[maxRow * width + _OFFSET],
								&u[_OFFSET * width + _OFFSET]);
							RowOperation<_T, _OFFSET>::swap(&l[maxRow * width],
								&l[_OFFSET * width]);
						}

						return decompose(u, l);
					}

					template<uint32_t _BXWIDTH>
					static bool solveLyEqualsB(Matrix<_HEIGHT, _WIDTH, _T>* l,
						Matrix<_HEIGHT, _BXWIDTH, _T>* bx)
					{
						_T factor = (*l)[_OFFSET][_OFFSET];
						RowOperation<_T, _BXWIDTH>::multiply((*bx)[_OFFSET],
							(*bx)[_OFFSET], 1.0 / factor);

						// make sure there is a row under us
						if (_OFFSET >= _HEIGHT - 1)
							return true;

						// substract the row so that all lower rows have a 0 at _OFFSET
						for (uint32_t j = _OFFSET + 1; j < _HEIGHT; ++j)
						{
							_T factor = (*l)[j][_OFFSET];
							RowOperation<_T, _BXWIDTH>::addRowTimesFactor(
								(*bx)[j], (*bx)[j], (*bx)[_OFFSET], -factor);
							(*l)[j][_OFFSET] = 0;
						}
						// solve the problem for _SIZE-1
						LUSubDecomposition<_T, _OFFSET + 1, _WIDTH,
							_HEIGHT>::solveLyEqualsB(l, bx);

						return true;
					}

					template<uint32_t _BXWIDTH>
					static bool solveUxEqualsY(Matrix<_HEIGHT, _WIDTH, _T>* u,
						Matrix<_HEIGHT, _BXWIDTH, _T>* bx)
					{
						// make sure there is a row under us
						if (_OFFSET >= _HEIGHT - 1)
							return true;

						// solve the problem for _SIZE-1
						LUSubDecomposition<_T, _OFFSET + 1, _WIDTH,
							_HEIGHT>::solveUxEqualsY(u, bx);

						// substract the row so that all upper rows have a 0 at _OFFSET
						for (uint32_t j = 0; j < _OFFSET + 1; ++j)
						{
							_T factor = (*u)[j][_OFFSET + 1];
							RowOperation<_T, _BXWIDTH>::addRowTimesFactor(
								(*bx)[j], (*bx)[j], (*bx)[_OFFSET + 1],
								-factor);
							(*u)[j][_OFFSET + 1] = 0;
						}

						return true;
					}

					template<uint32_t _BXWIDTH>
					static bool solve(const Matrix<_HEIGHT, _WIDTH, _T>& l,
						const Matrix<_HEIGHT, _WIDTH, _T>& u,
						Matrix<_HEIGHT, _BXWIDTH, _T>* bx)
					{
						// we want to solve Ax = b where A = LU
						// we can write LUx = b
						// subsituing Ux = y, we have:
						//				Ly = b
						//
						// Ux = y is easy to solve, and Ly = b is also easy to solve

						// start solving Ly = b

						Matrix<_HEIGHT, _WIDTH, _T> lCopy(l);
						if (!solveLyEqualsB(&lCopy, bx))
							return false;

						Matrix<_HEIGHT, _WIDTH, _T> uCopy(u);
						if (!solveUxEqualsY(&uCopy, bx))
							return false;

						return true;
					}
			};

			template<typename _T, uint32_t _SIZE>
			class RowOperation
			{
				public:
					static void
					multiply(_T* intoRow, const _T* row, const _T& factor)
					{
						for (uint32_t i = 0; i < _SIZE; ++i)
							intoRow[i] = row[i] * factor;
					}

					static void addRowTimesFactor(_T* intoRow, const _T* srcRow,
						const _T* addRow, const _T& timesFactor)
					{
						for (uint32_t i = 0; i < _SIZE; ++i)
							intoRow[i] = srcRow[i] + addRow[i] * timesFactor;
					}

					static void swap(_T* row1, _T* row2)
					{
						_T tmp[_SIZE];
						memcpy(tmp, row1, _SIZE * sizeof(_T));
						memcpy(row1, row2, _SIZE * sizeof(_T));
						memcpy(row2, tmp, _SIZE * sizeof(_T));
					}
			};

			template<typename _T>
			class RowOperation<_T, 0>
			{
				public:
					static void
					multiply(_T* intoRow, const _T* row, const _T& factor)
					{
					}

					static void addRowTimesFactor(_T* intoRow, const _T* srcRow,
						const _T* addRow, const _T& timesFactor)
					{
					}

					static void swap(_T* row1, _T* row2)
					{
					}
			};

			template<typename _T, uint32_t _OFFSET, uint32_t _WIDTH>
			class LUSubDecomposition<_T, _OFFSET, _WIDTH, _OFFSET>
			{
				public:
					static bool decomposeRecur(_T* u, _T* l, int* p)
					{
						return true;
					}

					static bool decomposeRecur(_T* u, _T* l)
					{
						return true;
					}

					template<uint32_t _BXWIDTH>
					static bool solveUxEqualsY(Matrix<_OFFSET, _WIDTH, _T>* u,
						Matrix<_OFFSET, _BXWIDTH, _T>* bx)
					{
						return true;
					}

					template<uint32_t _BXWIDTH>
					static bool solveLyEqualsB(Matrix<_OFFSET, _WIDTH, _T>* l,
						Matrix<_OFFSET, _BXWIDTH, _T>* bx)
					{
						return true;
					}
			};
	};
}

#endif
