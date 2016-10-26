cmake_minimum_required(VERSION 3.5)
project(geometry)
set(SOURCE_FILES GeometryCalculator.cpp GeometryCalculator.h Line.h
		LineSegment.h Ray.h SpatialTree.h Vector.h)
foreach (loop_var ${SOURCE_FILES})
	set(MY_SOURCE_FILES ${MY_SOURCE_FILES} ${PROJECT_NAME}/${loop_var})
endforeach(loop_var)
set(${PROJECT_NAME}_SOURCE_FILES ${MY_SOURCE_FILES} PARENT_SCOPE)