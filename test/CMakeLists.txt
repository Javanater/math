cmake_minimum_required(VERSION 3.5)
project(MathTests)
add_executable(${PROJECT_NAME}
		GeometryCalculatorTest.cpp
		ReferenceFrameTest.cpp
		MathTest.cpp)
target_link_libraries(${PROJECT_NAME} gtest gtest_main)
target_link_libraries(${PROJECT_NAME} Math)
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
