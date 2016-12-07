cmake_minimum_required(VERSION 3.5)
project(Math)
set(SOURCE_FILES
		Math.hpp

		#geometry
		geometry/Geometry.hpp
		geometry/GeometryCalculator.cpp
		geometry/GeometryCalculator.hpp
		geometry/Line.hpp
		geometry/LineSegment.hpp
		geometry/ReferenceFrame.hpp
		geometry/ReferencePoint.hpp
		geometry/Ray.hpp
		#		geometry/SpatialTree.h
		#[[geometry/Vector.h]]  )
add_library(${PROJECT_NAME} ${SOURCE_FILES})

