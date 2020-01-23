cmake_minimum_required(VERSION 3.5)

include_guard(GLOBAL)

project(Math)

set(SOURCE_FILES
    include/math/Math.hpp
    include/math/geometry/Geometry.hpp
    include/math/geometry/GeometryCalculator.cpp
    include/math/geometry/GeometryCalculator.hpp
    include/math/geometry/Line.hpp
    include/math/geometry/LineSegment.hpp
    include/math/geometry/ReferenceFrame.hpp
    include/math/geometry/ReferencePoint.hpp
    include/math/geometry/Ray.hpp
)

add_library(${PROJECT_NAME} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME}
	PUBLIC
		$<INSTALL_INTERFACE:include>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/src
)

find_package(Eigen3 REQUIRED)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)

find_package(Boost COMPONENTS unit_test_framework REQUIRED)
find_package(GTest REQUIRED)

add_executable(${PROJECT_NAME}_Test
	test/GeometryCalculatorTest.cpp
	test/ReferenceFrameTest.cpp
	test/MathTest.cpp
)

target_link_libraries(${PROJECT_NAME}_Test
	${PROJECT_NAME}
	Boost::unit_test_framework
    GTest::Main
)

