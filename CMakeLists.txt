cmake_minimum_required(VERSION 3.5)
project(Math)
add_subdirectory(geometry)
#add_subdirectory(matrix)
#set(SOURCE_FILES Math.hpp ${geometry_SOURCE_FILES} ${matrix_SOURCE_FILES} geometry/Geometry.hpp)
set(SOURCE_FILES Math.hpp ${geometry_SOURCE_FILES}
		#geometry
		geometry/Geometry.hpp
		#matrix
		matrix/matrix.cpp
		matrix/matrix.h
		matrix/LUDecomposition.h
		geometry/Pose.hpp)
add_library(${PROJECT_NAME} ${SOURCE_FILES})

