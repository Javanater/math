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
 * Helper function, not intended for external use
 *
 * @tparam T
 * @param t
 * @return
 */
template<class T>
inline T hypotH(const T t)
{
	return t * t;
}

/**
 * Helper function, not intended for external use
 *
 * @tparam T
 * @tparam Args
 * @param t
 * @param args
 * @return
 */
template<class T, class... Args>
inline T hypotH(const T t, const Args ... args)
{
	return t * t + hypotH<Args...>(args...);
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
	return std::sqrt(hypotH<T, Args...>(t, args...));
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
	double r = a - b;
	r = unsignedMod<T>(r + boost::math::constants::pi<T>(),
		boost::math::constants::two_pi<T>()) - boost::math::constants::pi<T>();
	return r;
}
}

#endif //PROJECTS_MATH_HPP
