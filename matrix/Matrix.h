/*
 * Matrix.h
 *
 *  Created on: Jun 20, 2016
 *      Author: Madison
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <initializer_list>
#include <array>
#include <ostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <utilities/utilities.h>
#include <utilities/Pretty.h>
#include <cmath>
#include "LUDecomposition.h"

namespace flabs
{
	template<uint32_t DIM, typename ValueType>
	class Vector;

	template<class Mat, uint32_t SKIP>
	struct Minor
	{
		typedef typename Mat::Type Type;
		static const uint32_t      rows   = Mat::rows;
		static const uint32_t      cols   = Mat::cols;
		static const uint32_t      ccount = Mat::ccount - 1;
		static const uint32_t      rstart = Mat::rstart + 1;
		const Mat& mat;

		Minor(const Mat& mat) : mat(mat)
		{
		}

		inline Type operator()(uint32_t row, uint32_t col) const
		{
			if (col >= SKIP)
				return mat(row + 1, col + 1);
			else
				return mat(row + 1, col);
		}
	};

	template<uint32_t i, class Mat>
	struct DetHelp;

	template<uint32_t i, class Mat>
	struct Det
	{
		static inline typename Mat::Type determinantBR(const Mat& mat)
		{
			return DetHelp<Mat::ccount - 1, Mat>::determinantBLoop(mat);
		}
	};

	template<class Mat>
	struct Det<1, Mat>
	{
		static inline typename Mat::Type determinantBR(const Mat& mat)
		{
			return mat(0, 0);
		}
	};

	template<class Mat>
	struct Det<0, Mat>
	{
		static inline typename Mat::Type determinantBR(const Mat& mat)
		{
			return 1;
		}
	};

	template<uint32_t i, class Mat>
	struct DetHelp
	{
		static inline typename Mat::Type determinantBLoop(const Mat& mat)
		{
			typename Mat::Type
				temp                =
				DetHelp<i - 1, Mat>::determinantBLoop(mat);

			typedef Minor<Mat, i> Mat2;
			Mat2                  minor(mat);
			typename Mat::Type    t = Det<Mat::rows - Mat2::rstart,
				Mat2>::determinantBR(minor);

			if (i % 2 == 0)
				return temp + mat(0, i) * t;
			else
				return temp - mat(0, i) * t;
		}
	};

	template<class Mat>
	struct DetHelp<0, Mat>
	{
		static inline typename Mat::Type determinantBLoop(const Mat& mat)
		{
			typedef Minor<Mat, 0> Mat2;
			Mat2                  minor(mat);
			typename Mat::Type    t = Det<Mat::rows - Mat2::rstart,
				Mat2>::determinantBR(minor);
			return mat(0, 0) * t;
		}
	};

	// n x m : ROWS x COLS
	template<uint32_t ROWS, uint32_t COLS, class TYPE = double>
	class Matrix
	{
		public:
			typedef Matrix<ROWS, COLS, TYPE> Mat;
			typedef TYPE                     Type;

			static const uint32_t rows   = ROWS;
			static const uint32_t cols   = COLS;
			static const uint32_t ccount = COLS;
			static const uint32_t rstart = 0;

			TYPE data[ROWS * COLS];

			Matrix()
			{
			}

			Matrix(const std::initializer_list<
				std::initializer_list<TYPE>>& initialValues)
			{
				uint32_t r = 0;
				for (const auto& a : initialValues)
				{
					uint32_t c = 0;
					for (const auto& b : a)
					{
						operator()(r, c) = b;
						++c;
					}
					++r;
				}
			}

		private:
			static inline void init(Mat& matrix)
			{
			}

			template<class Vector, class... Vectors>
			static inline void
			init(Mat& matrix, const Vector& vector, const Vectors& ... vectors)
			{
				for (uint32_t r = 0; r < ROWS; ++r)
					matrix(r, COLS - sizeof...(vectors) - 1) = vector[r];
				init(matrix, vectors...);
			}

		public:
			template<class... Vectors>
			Matrix(const Vector<ROWS, TYPE>& vector,
				const Vectors& ... vectors)
			{
				static_assert(sizeof...(vectors) == COLS - 1,
					"Must provide same number of vectors as columns.");
				init(*this, vector, vectors...);
			}

			Matrix(
				const std::array<std::array<TYPE, COLS>, ROWS>& initialValues)
			{
				for (uint32_t r = 0; r < ROWS; ++r)
					for (uint32_t c = 0; c < COLS; ++c)
						operator()(r, c) = initialValues[r][c];
			}

			Matrix(std::function<TYPE()>& generator)
			{
				for (uint32_t r = 0; r < ROWS; ++r)
					for (uint32_t c = 0; c < COLS; ++c)
						operator()(r, c) = generator();
			}

			Matrix(TYPE fill)
			{
				for (uint32_t r = 0; r < ROWS; ++r)
					for (uint32_t c = 0; c < COLS; ++c)
						operator()(r, c) = fill;
			}

			~Matrix()
			{
			}

			static Mat identity()
			{
				Matrix<ROWS, COLS, TYPE> result;
				for (uint32_t            r = 0; r < ROWS; ++r)
					for (uint32_t c = 0; c < COLS; ++c)
						if (c == r)
							result(r, c) = 1;
						else
							result(r, c) = 0;
				return result;
			}

			template<uint32_t MCOLS>
			inline Matrix<ROWS, MCOLS, TYPE>
			operator*(const Matrix<COLS, MCOLS, TYPE>& matrix) const
			{
				Matrix<ROWS, MCOLS, TYPE> result;
				for (uint32_t             r = 0; r < ROWS; ++r)
					for (uint32_t c = 0; c < MCOLS; ++c)
					{
						result(r, c) = 0;
						for (uint32_t i = 0; i < COLS; ++i)
							result(r, c) += operator()(r, i) * matrix(i, c);
					}
				return result;
			}

			inline Type& operator()(uint32_t row, uint32_t col)
			{
				return data[row * COLS + col];
			}

			inline const Type& operator()(uint32_t row, uint32_t col) const
			{
				return data[row * COLS + col];
			}

			Type* operator[](int row)
			{
				return &data[row * COLS];
			}

			const Type* operator[](int row) const
			{
				return &data[row * COLS];
			}

			inline TYPE det() const
			{
				return DetHelp<ccount - 1, Mat>::determinantBLoop(*this);
			}

			inline Mat inv() const
			{
				Mat l, u;
				LUDecomposition::decompose(*this, &l, &u);
				Mat i = identity();
				LUDecomposition::solve(l, u, &i);
				return i;
			}

			//TODO: Increase efficiency
			inline bool equals(const Mat& matrix,
				double percentageTolerance = .00000001) const
			{
				for (uint32_t r = 0; r < ROWS; ++r)
					for (uint32_t c = 0; c < COLS; ++c)
						if (std::abs(operator()(r, c) - matrix(r, c)) /
							operator()(r, c) > percentageTolerance)
							return false;
				return true;
			}

			inline Mat& operator*=(const TYPE scaler)
			{

				for (uint32_t r = 0; r < ROWS; ++r)
					for (uint32_t c = 0; c < COLS; ++c)
						operator()(r, c) *= scaler;
				return *this;
			}

			inline Mat& operator/=(const TYPE scaler)
			{

				for (uint32_t r = 0; r < ROWS; ++r)
					for (uint32_t c = 0; c < COLS; ++c)
						operator()(r, c) /= scaler;
				return *this;
			}

			inline Mat& operator+=(const TYPE scaler)
			{

				for (uint32_t r = 0; r < ROWS; ++r)
					for (uint32_t c = 0; c < COLS; ++c)
						operator()(r, c) += scaler;
				return *this;
			}

			inline Mat& operator-=(const TYPE scaler)
			{

				for (uint32_t r = 0; r < ROWS; ++r)
					for (uint32_t c = 0; c < COLS; ++c)
						operator()(r, c) -= scaler;
				return *this;
			}

			inline Matrix<COLS, ROWS, TYPE> transposed() const
			{
				Matrix<COLS, ROWS, TYPE> matrix;
				for (uint32_t            r = 0; r < COLS; ++r)
					for (uint32_t c = 0; c < ROWS; ++c)
						matrix(r, c) = operator()(c, r);
				return matrix;
			}

			inline Mat& transpose()
			{
				return *this = transposed();
			}

			inline TYPE* ptr()
			{
				return (TYPE*) data;
			}

			inline const TYPE* ptr() const
			{
				return (TYPE*) data;
			}

			inline bool isUpper() const
			{
				for (uint32_t r = 1; r < ROWS; ++r)
					for (uint32_t c = 0; c < r; ++c)
						if (operator()(r, c) != 0)
							return false;
				return true;
			}

			inline bool isLower() const
			{
				for (uint32_t r = 0; r < ROWS; ++r)
					for (uint32_t c = r + 1; c < COLS; ++c)
						if (operator()(r, c) != 0)
							return false;
				return true;
			}

			inline bool isSingular(
				TYPE tolerance = std::numeric_limits<TYPE>::epsilon() * 4) const
			{
				return det() <= tolerance;
			}

			inline Mat& operator=(std::function<TYPE()>& generator)
			{
				for (uint32_t r = 0; r < ROWS; ++r)
					for (uint32_t c = 0; c < COLS; ++c)
						operator()(r, c) = generator();
				return *this;
			}

			template<uint32_t R1, uint32_t C1,
				class T1, uint32_t R2, uint32_t C2, class T2>
			friend inline bool
			operator==(Matrix<R1, C1, T1> matrix1, Matrix<R2, C2, T2> matrix2);

			template<uint32_t R1, uint32_t C1, class T1>
			friend std::ostream&
			operator<<(std::ostream& out, const Matrix<R1, C1, T1>& mat);
	};

	template<uint32_t R1, uint32_t C1, class T1, uint32_t R2, uint32_t C2,
		class T2>
	inline bool
	operator==(Matrix<R1, C1, T1> matrix1, Matrix<R2, C2, T2> matrix2)
	{
		if (R1 != R2 || C1 != C2)
			return false;

		for (uint32_t r = 0; r < R1; ++r)
			for (uint32_t c = 0; c < C1; ++c)
				if (matrix1(r, c) != matrix2(r, c))
					return false;

		return true;
	}

	template<uint32_t R1, uint32_t C1, class T1>
	std::ostream& operator<<(std::ostream& out, const Matrix<R1, C1, T1>& mat)
	{
		out << "{";

		if (mat.rows > 0 && mat.cols > 0)
		{
			out << "{" << mat(0, 0);

			for (uint32_t c = 1; c < mat.cols; ++c)
				out << ", " << mat(0, c);

			out << "}";

			for (uint32_t r = 1; r < mat.rows; ++r)
			{
				out << ", {" << mat(r, 0);

				for (uint32_t c = 1; c < mat.cols; ++c)
					out << ", " << mat(r, c);

				out << "}";
			}
		}

		out << "}";

		return out;
	}

	template<uint32_t ROWS, uint32_t COLS>
	class Pretty<Matrix<ROWS, COLS>>
	{
		private:
			typedef Matrix<ROWS, COLS> Type;
			const Type& t;

		public:
			Pretty(const Type& t) : t(t)
			{
			}

			~Pretty()
			{
			}

			friend std::ostream&
			operator<<(std::ostream& out, const Pretty<Type>& pretty)
			{
				constexpr uint32_t rows = Type::rows - Type::rstart;
				constexpr uint32_t cols = Type::ccount;
				std::string        stringMatrix[rows][cols];

				for (uint32_t r = 0; r < rows; ++r)
					for (uint32_t c = 0; c < cols; ++c)
						stringMatrix[r][c] = toString(pretty.t(r, c));

				uint32_t      columnWidths[cols];
				for (uint32_t c            = 0; c < cols; ++c)
					columnWidths[c] = 0;

				for (uint32_t r = 0; r < rows; ++r)
					for (uint32_t c = 0; c < cols; ++c)
						columnWidths[c] = std::max(columnWidths[c],
							(uint32_t) stringMatrix[r][c].length());

				uint32_t      totalWidth = 0;
				for (uint32_t c          = 0; c < cols; ++c)
					totalWidth += columnWidths[c];

				for (uint32_t r = 0; r < rows; ++r)
				{
					out << "[";
					out << stringMatrix[r][0];
					for (uint32_t c = 1; c < cols; ++c)
					{
						out << ' ';
						out.width(columnWidths[c]);
						out << stringMatrix[r][c];
						out.width(0);
					}
					out << "]\n";
				}

				return out;
			}
	};
}

#endif
