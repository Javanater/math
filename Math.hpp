//
// Created by Madison on 8/24/2016.
//

#ifndef PROJECTS_MATH_HPP
#define PROJECTS_MATH_HPP

#include "geometry/Geometry.hpp"
#include <boost/math/constants/constants.hpp>

namespace flabs
{
/**
 * Computes a^2
 *
 * @tparam T
 * @param t
 * @return
 */
template<class T>
inline T sumSquares(const T t)
{
	return t * t;
}

/**
 * Computes a^2 + b^2 + ...
 *
 * @tparam T
 * @tparam Args
 * @param t
 * @param args
 * @return
 */
template<class T, class... Args>
inline T sumSquares(const T t, const Args ... args)
{
	return t * t + sumSquares<Args...>(args...);
}

/**
 * Computes sqrt(a^2 + b^2 + ...)
 *
 * @tparam T
 * @tparam Args
 * @param t
 * @param args
 * @return
 */
template<class T, class... Args>
inline T hypot(const T t, const Args ... args)
{
	return std::sqrt(sumSquares<T, Args...>(t, args...));
}

/**
 * Computes the unsigned floating point remainder of a / b
 *
 * @tparam T
 * @param a
 * @param b
 * @return
 */
template<class T>
inline T unsignedMod(T a, T b)
{
	return a - floor(a / b) * b;
}

template<class T>
inline T zeroTo2Pi(T a)
{
	return unsignedMod<T>(a, boost::math::constants::two_pi<T>());
}

/**
 * Computes the smallest difference between angles a and b within the range
 * [-pi, pi).
 *
 * @tparam T: the data type
 * @param a: first angle in radians
 * @param b: second angle in radians
 * @return a - b [-pi, pi)
 */
template<class T>
inline T angleDifference(T a, T b)
{
	a   = zeroTo2Pi<T>(a);
	b   = zeroTo2Pi<T>(b);
	T r = a - b;
	if (r >= boost::math::constants::pi<T>())
		r -= boost::math::constants::two_pi<T>();
	else if (r < -boost::math::constants::pi<T>())
		r += boost::math::constants::two_pi<T>();
	return r;
}

/**
 * Computes the smallest difference between inputs a and b within the range
 * [min, max).
 *
 * @tparam T: the data type
 * @tparam min: minimum of interval inclusive
 * @tparam max: maximum of interval exclusive
 * @param a: first measure
 * @param b: second measure
 * @return a - b [min, max)
 */
template<class T>
inline T intervalDifference(T a, T b, T min, T max)
{
	T r = a - b;
	r = unsignedMod<T>(r - min, max - min) + min;
	return r;
}

/**
 * Computes the smallest difference between inputs a and b within the range
 * [min, max).
 *
 * @tparam T: the data type
 * @tparam min: minimum of interval inclusive
 * @tparam max: maximum of interval exclusive
 * @param a: first measure
 * @param b: second measure
 * @return a - b [min, max)
 */
template<class T, T min, T max>
inline T intervalDifference(T a, T b)
{
	T r = a - b;
	r = unsignedMod<T>(r - min, max - min) + min;
	return r;
}
}

#endif //PROJECTS_MATH_HPP
